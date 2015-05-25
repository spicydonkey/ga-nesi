#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include "AdvXMLParser.h"
#include "GAEngine.h"
#include "cellml-api-cxx-support.hpp"
#include "IfaceCellML_APISPEC.hxx"
#include "IfaceCIS.hxx"
#include "CellMLBootstrap.hpp"
#include "CISBootstrap.hpp"
#include "virtexp.h"
#include "distributor.h"


using namespace std;
using namespace AdvXMLParser;
using namespace iface::cellml_api;
using namespace iface::cellml_services;

//If SUPPORT_BLOCK_SAMPLING is defined
//the selection algorithms may use either
//probabilistic search or blocking search
//#define SUPPORT_BLOCK_SAMPLING

ObjRef<iface::cellml_api::CellMLBootstrap> bootstrap; //CellML api bootstrap
ObjRef<iface::cellml_services::CellMLIntegrationService> cis;

VariablesHolder var_template; //template for the variables, just holds names of the variables

int verbosity=0;

void usage(const char *name)
{
    printf("Usage: %s <experiment definition xml> [-v [-v]]\n",name);
    printf("Where -v increases the verbosity of the output\n");
}

//Open and read XML configuration file
char *OpenXmlFile(const char *name,long& nSize)
{
    FILE *f=fopen(name,"rb");
    char *pBuffer=NULL;

	//check for file open error
    if(!f)
        return NULL;
	
	//obtain file size
    fseek(f,0,SEEK_END);
    nSize=ftell(f);
    fseek(f,0,SEEK_SET);

	//allocate memory to contain the whole file
    pBuffer=new char[nSize+1];

	//copy the file into buffer
    fread(pBuffer,nSize,1,f);	// fread usage!?
	//fread(pBuffer,1,nSize,f);
    pBuffer[nSize]=0;
    fclose(f);

    return pBuffer;
}

//Initialise GA engine
int SetAndInitEngine(GAEngine<COMP_FUNC >& ga,const AdvXMLParser::Element& elem)
{
	//Get GA parameters from XML file
    int initPopulation=atoi(elem.GetAttribute("InitialPopulation").GetValue().c_str());
    double mutation=atof(elem.GetAttribute("Mutation_proportion").GetValue().c_str());
    double cross=atof(elem.GetAttribute("Crossover_proportion").GetValue().c_str());
    int generations=atoi(elem.GetAttribute("Generations").GetValue().c_str());
    int block_sample=atoi(elem.GetAttribute("Sampling").GetValue().c_str());
    

    //Set defaults for uninitialised variables
    if(!initPopulation)
        initPopulation=100;

    if(cross>1.0)
        cross=1.0;
    if(mutation>1.0)
        mutation=1.0;
    ga.prob_cross()=cross;
    ga.prob_mutate()=mutation;
    ga.part_cross()=(int)((double)initPopulation*cross);
    ga.part_mutate()=(int)((double)initPopulation*mutation);
#ifdef SUPPORT_BLOCK_SAMPLING
    ga.block_sample()=(block_sample==0);
#endif
    
    //Read alleles information
    for(int i=0;;i++)
    {
        const AdvXMLParser::Element& al=elem("Alleles",0)("Allele",i);		// what does this line do??
        std::wstring name; 
        if(al.IsNull())
           break;
        name=convert(al.GetAttribute("Name").GetValue());
        ga.AddAllele(name);
        //Set allele limits
        ga.AddLimit(name,atof(al.GetAttribute("LowerBound").GetValue().c_str()),atof(al.GetAttribute("UpperBound").GetValue().c_str()));
        //Initialise variable template
        var_template(name,0.0);
    }
    ga.set_borders(initPopulation);
    return (generations?generations:1);
}

void initialize_template_var(const AdvXMLParser::Element& elem)
{
    for(int i=0;;i++)
    {
        const AdvXMLParser::Element& al=elem("Alleles",0)("Allele",i);
        std::wstring name; 
        if(al.IsNull())
           break;
        name=convert(al.GetAttribute("Name").GetValue());
        var_template(name,0.0);
    }
}

//Observer callback
bool observer(WorkItem *w,double answer,void *g)
{
    GAEngine<COMP_FUNC> *ga=(GAEngine<COMP_FUNC> *)g;
   
    ga->process_workitem(w,answer);
    return true;
}



double do_compute(std::vector<double>& val)
{
    var_template.fillup(val);
    return VEGroup::instance().Evaluate(var_template);
}

//Slave process
//Returns only when quit command is received from the master
void run_slave(int proc)
{
    double req;
    MPI_Status stat;
    std::vector<double> data;

    var_template.collate(data); 
    while(1)
    {
        //check if data is received
        MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        if(stat.MPI_TAG==TAG_QUIT)
        {
            //Quit signal received
            break;
        }
        //Receive compute request and process it
        MPI_Recv(&data[0],data.size(),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        req=do_compute(data);
        //returns the result of the computations
        MPI_Send(&req,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    }
}

int main(int argc,char *argv[])
{
    char *pBuffer=NULL;
    long nSize=0;
    GAEngine<COMP_FUNC > ga;
    int proc,nproc;
    int generations=1;
    const char *filename=NULL;

    srand(time(NULL));
    MPI_Init(&argc,&argv);
    if(argc<2)
    {
        usage(argv[0]);
        return -1;
    }

    for(int i=1;i<argc;i++)
    {
        if(!strcmp(argv[i],"-v"))
            verbosity++;
        else
            filename=argv[i];
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    //Load and initialise CellML API
    bootstrap=CreateCellMLBootstrap();
    cis=CreateIntegrationService();

    if((pBuffer=OpenXmlFile(filename,nSize))==NULL)
    {
         fprintf(stderr,"Error opening input file %s\n",argv[1]);
         return -1;
    }

    //Load the experiments package
    try
    {
        Parser parser;
        ObjRef<CellMLComponentSet> comps;


        auto_ptr<Document> pDoc(parser.Parse(pBuffer,nSize));
        const Element& root=pDoc->GetRoot();
        for(int i=0;;i++)
        {
            VariablesHolder params;

            VirtualExperiment *vx=VirtualExperiment::LoadExperiment(root("VirtualExperiments",0)("VirtualExperiment",i));
            if(!vx)
               break;
            VEGroup::instance().add(vx);
        }
        if(!proc)
        {
            generations=SetAndInitEngine(ga,root("GA",0));
        }
        else
        {
            initialize_template_var(root("GA",0));
        }
    }
    catch(ParsingException e)
    {
        printf("Parsing error at line %d\n",e.GetLine());
    }
    delete [] pBuffer;
    //Wait until all the clints are ready
    //
    MPI_Barrier(MPI_COMM_WORLD);

    //Only master tasks needs GA engine to be initialised and used   
    if(!proc)
    {
        //Master task
        VariablesHolder v;


        ga.Initialise();
        ga.RunGenerations(generations); //Run generations
        double bf=ga.GetBest(v);
        printf("Best fitness: %lf\n",bf);
        //Print out results for best fitness
        for(int i=0;;i++)
        {
            wstring name=v.name(i);
            if(!name.size())
                break;
            printf("Best[%s]=%lf\n",convert(name).c_str(),v(name));
        }
        Distributor::instance().finish();
    }
    else
    {
        run_slave(proc);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}

