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

int verbosity=0;	// verbosity initialised to 0

void usage(const char *name)
{
    printf("Usage: %s <experiment definition xml> [-v [-v]]\n",name);
    printf("Where -v increases the verbosity of the output\n");
}

//Open and read XML configuration file
//nSize is assigned the size of file
//return the contents of the file as a null-terminated char array
char *OpenXmlFile(const char *name,long& nSize)
{
    FILE *f=fopen(name,"rb");	// open file "name" for reading
    char *pBuffer=NULL;	// initialise a buffer for storing C-string to a nullptr

	//check for file open error
    if(!f)
        return NULL;
	
	//obtain file size
    fseek(f,0,SEEK_END);
    nSize=ftell(f);
    fseek(f,0,SEEK_SET);

	//allocate memory to contain the whole file
    pBuffer=new char[nSize+1];

	
    fread(pBuffer,nSize,1,f);	//copy the file into buffer (usage of fread)
	//fread(pBuffer,1,nSize,f);
    
	pBuffer[nSize]=0;	// null terminate the char array buffer
    fclose(f);

    return pBuffer;
}

//Initialise GA engine
	//return number of generations to run GA
int SetAndInitEngine(GAEngine<COMP_FUNC >& ga,const AdvXMLParser::Element& elem)
{
	//Get GA parameters from XML file
    int initPopulation=atoi(elem.GetAttribute("InitialPopulation").GetValue().c_str());
    double mutation=atof(elem.GetAttribute("Mutation_proportion").GetValue().c_str());
    double cross=atof(elem.GetAttribute("Crossover_proportion").GetValue().c_str());
    int generations=atoi(elem.GetAttribute("Generations").GetValue().c_str());
    int block_sample=atoi(elem.GetAttribute("Sampling").GetValue().c_str());
    

    //Set the parameters for the GA engine accordingly
    if(!initPopulation)
        initPopulation=100;		// default population size if in XML it is undeclared or 0
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
        const AdvXMLParser::Element& al=elem("Alleles",0)("Allele",i);		// assign to "al" the ith Allele sub-element of the (first) Alleles element in "elem"
        std::wstring name;
        if(al.IsNull())
           break;	// exhaustively iterate through the "Allele" sub-elements of Alleles
        name=convert(al.GetAttribute("Name").GetValue());	// get the String value for Name attribute in "al" allele sub-element; convert it to wstring and store in name
        // Add allele
		ga.AddAllele(name);	
        // Set allele limits
        ga.AddLimit(name,atof(al.GetAttribute("LowerBound").GetValue().c_str()),atof(al.GetAttribute("UpperBound").GetValue().c_str()));
        // Initialise variable template
        var_template(name,0.0);
    }
    ga.set_borders(initPopulation);		// set max population of GA

	// return the num of generations to run GA
    return (generations?generations:1);	// default number of generations to run is 1
}


//Initialise the template variable holder with Alleles sub-element
void initialize_template_var(const AdvXMLParser::Element& elem)
{
    for(int i=0;;i++)
    {
        const AdvXMLParser::Element& al=elem("Alleles",0)("Allele",i);	// assign to "al" the ith Allele sub-element of the (first) Alleles element in "elem"
        std::wstring name; 
        if(al.IsNull())
           break;	// exhaustively iterate through the "Allele" sub-elements of Alleles
        name=convert(al.GetAttribute("Name").GetValue());	// get the String value for Name attribute in "al" allele sub-element; convert it to wstring and store in name
        // Initialise variable template
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


// perform Evaluate from given vector of doubles
double do_compute(std::vector<double>& val)
{
	// fill-up the tmp's allele values with supplied data
    var_template.fillup(val);
	// evaluate this chromosome's fit and return the representative residual
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
    GAEngine<COMP_FUNC > ga;	// initialise a GA engine
    int proc,nproc;
    int generations=1;
    const char *filename=NULL;

    srand(time(NULL));	// seed the RNG

    MPI_Init(&argc,&argv);

    if(argc<2)
    {
		// warn user for incorrect usage of command
        usage(argv[0]);
        return -1;
    }

    for(int i=1;i<argc;i++)
    {
        if(!strcmp(argv[i],"-v"))
			// if an arg string is "-v" increment verbosity
            verbosity++;
        else
			// other arg string becomes the filename
            filename=argv[i];
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    //Load and initialise CellML API
    bootstrap=CreateCellMLBootstrap();
    cis=CreateIntegrationService();

	// Read input file and store contents in buffer
    if((pBuffer=OpenXmlFile(filename,nSize)) == NULL)
    {
		fprintf(stderr,"Error opening input file %s\n",argv[1]);
		return -1;
    }
	// Read success: pBuffer is a C-string containing file and nSize is size of file

    //Load the experiments package: GA parameters and Virtual Experiment data
    try
    {
        Parser parser;
        ObjRef<CellMLComponentSet> comps;	//??? comps not referenced elsewhere in project

		// Parse the XML contents in buffer
        auto_ptr<Document> pDoc(parser.Parse(pBuffer,nSize));	// can throw an exception
		// Get the root of the XML structure
        const Element& root=pDoc->GetRoot();

		
		// load all virtual experiments in the XML file
        //
		for(int i=0;;i++)
        {
            VariablesHolder params;	//??? unused

			// load the ith VE in file
            VirtualExperiment *vx=VirtualExperiment::LoadExperiment(root("VirtualExperiments",0)("VirtualExperiment",i));
			// check if this VE is defined
			if(!vx)
               break;
			
			// add each VE into the group singleton
            VEGroup::instance().add(vx);
        }

		
		// load the GA parameters from file and initialise the engine
		//
        if(!proc)
        {
			// assign number of generations and initialise the parameters for the GA engine
            generations=SetAndInitEngine(ga,root("GA",0));
        }
        else
        {
			// initialise the template variable holder
            initialize_template_var(root("GA",0));
        }
    }
    catch(ParsingException e)
    {
        printf("Parsing error at line %d\n",e.GetLine());
    }
    delete [] pBuffer;	// free memory used to store file

	//Wait until all the clints are ready
    //
    MPI_Barrier(MPI_COMM_WORLD);

    //Only master tasks needs GA engine to be initialised and used   
    if(!proc)
    {
        //Master task
        VariablesHolder v;

		//Initialise the population in GA engine
        ga.Initialise();
		//Run GA
        ga.RunGenerations(generations);
        
		double bf=ga.GetBest(v);	// v stores the best Genome's chromosome from the run; bf stores its fitness
        
		//Print out results for best fitness
		printf("Best fitness: %lf\n",bf);
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

