#include <unistd.h>
#include <stdlib.h>
#include "virtexp.h"
#include "AdvXMLParser.h"
#include "utils.h"
#include "cellml_observer.h"
#include <math.h>


using namespace std;
using namespace AdvXMLParser;
using namespace iface::cellml_api;
using namespace iface::cellml_services;


#define EPSILON 0.01

extern ObjRef<iface::cellml_api::CellMLBootstrap> bootstrap; //CellML api bootstrap
extern ObjRef<iface::cellml_services::CellMLIntegrationService> cis;

VirtualExperiment::VirtualExperiment():m_nResultColumn(-1),m_ReportStep(0.0),m_MaxTime(0),m_Accuracy(EPSILON)
{
}

VirtualExperiment::~VirtualExperiment()
{
}

VirtualExperiment *VirtualExperiment::LoadExperiment(const AdvXMLParser::Element& elem)
{
    VirtualExperiment *vx=NULL;

    string strName=elem.GetAttribute("ModelFilePath").GetValue();
	
	// Check if the model file path is defined
    if(!strName.size())
        return NULL;
    vx=new VirtualExperiment;
    if(!vx->LoadModel(strName))
    {
        delete vx;
        vx=NULL;
    }
    else //read assessment points
    {
        vx->m_nResultColumn=atoi(elem.GetAttribute("ResultColumn").GetValue().c_str());

        if(elem.GetAttribute("Accuracy").GetValue().size())
              vx->m_Accuracy=atof(elem.GetAttribute("Accuracy").GetValue().c_str());
       
        vx->m_MaxTime=atoi(elem.GetAttribute("MaxSecondsForSimulation").GetValue().c_str());
        vx->m_ReportStep=atof(elem.GetAttribute("ReportStep").GetValue().c_str());
        for(int i=0;;i++)
        {
            const AdvXMLParser::Element& al=elem("AssessmentPoints",0)("AssessmentPoint",i);
            POINT pt;

            if(al.IsNull())
                break;
            pt.first=atof(al.GetAttribute("time").GetValue().c_str());
            pt.second=atof(al.GetAttribute("target").GetValue().c_str());
            vx->m_Timepoints.push_back(pt);
        }
        //read parameters
        for(int i=0;;i++)
        {
            const AdvXMLParser::Element& al=elem("Parameters",0)("Parameter",i);

            if(al.IsNull())
                break;
            wstring name=convert(al.GetAttribute("ToSet").GetValue());
            double val=atof(al.GetAttribute("Value").GetValue().c_str());
            if(name.size())
                vx->m_Parameters[name]=val;
        }
    }
    
    return vx;
}

double VirtualExperiment::getSSRD(std::vector<std::pair<int,double> >& d)
{
    double r=0.0;

    for(int i=0;i<d.size();i++)
    {
        r+=pow((d[i].second-m_Timepoints[d[i].first].second)/m_Timepoints[d[i].first].second,2);
    }
    return r;
}

bool VirtualExperiment::LoadModel(const std::string& model_name)
{
    bool res=false;

    std::wstring modelURL=convert(model_name);
    m_strModelName=model_name;
 
    try
    {
        m_Model=bootstrap->modelLoader()->loadFromURL(modelURL); 
        res=true;
    }
    catch(CellMLException e)
    {
        printf("Error loading model %s\n",model_name.c_str());
        res=false;
    }
    return res;
}

void VirtualExperiment::SetParameters(VariablesHolder& v)
{
    for(int i=0;;i++)
    {
        wstring n=v.name(i);
        if(n.empty())
           break;
        double val=v(n);
        m_Parameters[n]=val; 
    }
}

void VirtualExperiment::SetVariables(VariablesHolder& v)
{
    ObjRef<iface::cellml_api::CellMLComponentSet> comps=m_Model->modelComponents();
    ObjRef<iface::cellml_api::CellMLComponentIterator> comps_it=comps->iterateComponents();
    ObjRef<iface::cellml_api::CellMLComponent> firstComp=comps_it->nextComponent();



    while(firstComp)
    {
        ObjRef<iface::cellml_api::CellMLVariableSet> vars=firstComp->variables();
        ObjRef<iface::cellml_api::CellMLVariableIterator> vars_it=vars->iterateVariables();
        ObjRef<iface::cellml_api::CellMLVariable> var=vars_it->nextVariable();


        string compname=convert(firstComp->name());

        while(var)
        {
            wstring name=var->name();
            wstring fullname=name;
            if(compname!="all" && compname!="")
            {
                fullname=convert(compname)+convert(".");
                fullname+=name;
            }
            if(v.exists(fullname))
            {
                char sss[120];
                gcvt(v(fullname),25,sss);
                std::wstring wv=convert(sss);
                var->initialValue(wv);
            }
            else if(m_Parameters.find(fullname)!=m_Parameters.end())
            {
                char sss[120];
                gcvt(m_Parameters[fullname],25,sss);
                std::wstring wv=convert(sss);
                var->initialValue(wv);
            }
            var=vars_it->nextVariable();
        }
        firstComp=comps_it->nextComponent();
    }
}

double VirtualExperiment::Evaluate()
{
    double res=0.0;
    //int j=0;
    ObjRef<iface::cellml_services::ODESolverCompiledModel> compiledModel;
    ObjRef<iface::cellml_services::ODESolverRun> osr;
    time_t calc_started;

    try
    {
       compiledModel=cis->compileModelODE(m_Model);
       osr=cis->createODEIntegrationRun(compiledModel);
       LocalProgressObserver *po=new LocalProgressObserver(compiledModel);

       osr->setProgressObserver(po);
       po->release_ref();
       osr->stepType(iface::cellml_services::BDF_IMPLICIT_1_5_SOLVE);
       osr->setStepSizeControl(1e-6,1e-6,1.0,0.0,1.0);
       osr->setResultRange(0.0,m_Timepoints[m_Timepoints.size()-1].first,m_Timepoints[m_Timepoints.size()-1].first);
       if(m_ReportStep)
            osr->setTabulationStepControl(m_ReportStep,true);

       calc_started=time(NULL);
       osr->start();
       while(!po->finished())
       {
           usleep(1000);
           if(m_MaxTime)
           {
               time_t t=time(NULL);
               if((unsigned long)(t-calc_started)>m_MaxTime)
               {
                   po->failed("Took too long to integrate");
                   throw CellMLException();
               }
           }
       }      

       if(!po->failed())
       {
           std::vector<double> vd;
           std::vector<std::pair<int,double> > results;
           int recsize=po->GetResults(vd);
 
           for(int i=0;i<vd.size();i+=recsize)
           {
               for(int j=0;j<m_Timepoints.size();j++)
                   if(in_range(vd[i],m_Timepoints[j].first,EPSILON))
                   {
                       results.push_back(make_pair(j,vd[i+m_nResultColumn]));
                       break;
                   }
           }
           if(!results.size())
           {
               fprintf(stderr,"Results vector is empty, Observer returned %d bytes (%d records)\n",vd.size(),vd.size()/recsize);
           }
           res=(results.size()?getSSRD(results):INFINITY);
       }
       else
           res=INFINITY;
    }
    catch(CellMLException e)
    {
        //fprintf(stderr,"Error evaluating model\n");
        res=INFINITY;
    }
    return res;
}


double VirtualExperiment::Runner::operator()(VariablesHolder& v)
{
    pOwner->SetVariables(v);
    return pOwner->Evaluate();
}


VEGroup::VEGroup()
{
}

VEGroup::~VEGroup()
{
}


VEGroup& VEGroup::instance()
{
    static VEGroup *pInstance=new VEGroup;

    return *pInstance;
}


double VEGroup::Evaluate(VariablesHolder& v)
{
    double res=0.0;
    int count=0;

    if(!experiments.size())
        return 0.0;
    for(int i=0;i<experiments.size();i++)
    {
        experiments[i]->SetVariables(v);
        double d=experiments[i]->Evaluate();
        if(d!=INFINITY)
        {
            res+=d;
            count++;
        }
    }
    return (count==experiments.size()?res/(double)count:INFINITY);
}


void VEGroup::add(VirtualExperiment *p)
{
    experiments.push_back(p);
}

