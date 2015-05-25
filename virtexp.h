#ifndef VIRTUAL_EXPERIMENT_H
#define VIRTUAL_EXPERIMENT_H
#include <stddef.h>
#include "cellml-api-cxx-support.hpp"
#include "IfaceCellML_APISPEC.hxx"
#include "CellMLBootstrap.hpp"
#include "AdvXMLParser.h"
#include "CISBootstrap.hpp"
#include "utils.h"
#include <string>
#include <functional>
#include <algorithm>

#define COMP_FUNC std::less_equal<double>

// define ALLELE as vector of <wstring, double> pairs
typedef std::vector<std::pair<std::wstring,double> > ALLELE;


/**
 *	VariablesHolder
 *	
 *	+ m_Vars : ALLELE
 *	
 *	
 **/
class VariablesHolder
{
    public:
	VariablesHolder() {}
	VariablesHolder(const VariablesHolder& other) { m_Vars.assign(other.m_Vars.begin(),other.m_Vars.end()); }	//copy other VarHolder's variables into this one
	~VariablesHolder() {}

	// = assign other VarHold to this if different 
	VariablesHolder& operator=(const VariablesHolder& other)
	{
	    if(&other!=this)
	       m_Vars.assign(other.m_Vars.begin(),other.m_Vars.end());
	    return *this;
	}

	//indexing VarHold obj by name:value in m_Vars
	double operator()(const std::wstring& name)
	{
	    ALLELE::iterator it=find_if(m_Vars.begin(),m_Vars.end(),
		      bind1st(pair_equal_to<std::wstring,double>(),name));	//find iterator to the pair in m_Vars for which the "first" member equals name (end if no such pair)
	    return (it==m_Vars.end()?double(0.0):it->second);
	}

	//Update an allele (pair) in VarHold's m_Var obj and return updated allele value
	double operator()(const std::wstring& name,double val)
	{
		//find the iterator for name
	    ALLELE::iterator it=find_if(m_Vars.begin(),m_Vars.end(),
		   bind1st(pair_equal_to<std::wstring,double>(),name));

	    if(it!=m_Vars.end())
	       it->second=val;	//update the allele if it already exists
	    else
	       m_Vars.push_back(std::make_pair<std::wstring,double>(std::wstring(name),double(val)));	//add the allele into m_Var if not stored yet
	    return val;
	}

	// TODO

	std::wstring name(int index)
	{
	    return ((index>=0 && index<m_Vars.size())?m_Vars[index].first:std::wstring());
	}

	bool exists(const std::wstring& name)
	{
	    ALLELE::iterator it=find_if(m_Vars.begin(),m_Vars.end(),
		   bind1st(pair_equal_to<std::wstring,double>(),name));
	    return (it!=m_Vars.end());
	}
        size_t size()
        {
            return m_Vars.size();
        }
	void collate(std::vector<double>& v)
	{
	     for(ALLELE::iterator it=m_Vars.begin();it!=m_Vars.end();++it)
	     {
			v.push_back(it->second);
	     }
	}
	void print()
	{
	     for(ALLELE::iterator it=m_Vars.begin();it!=m_Vars.end();++it)
	     {
		 printf("%s->%lf\n",convert(it->first).c_str(),it->second);
	     }
	}


	bool fillup(std::vector<double>& v)
	{
             if(v.size()!=m_Vars.size())
                 return false;
	     for(int i=0;i<v.size();i++)
	     {
		 m_Vars[i].second=v[i];
	     }
             return true;
	}

    private:
	ALLELE m_Vars;
};

class VirtualExperiment
{
    public:
        VirtualExperiment();
        ~VirtualExperiment();
        bool LoadModel(const std::string& model_name);
        static VirtualExperiment *LoadExperiment(const AdvXMLParser::Element& elem);
        void SetVariables(VariablesHolder& v);
        void SetParameters(VariablesHolder& v);
        double Evaluate();

        int resultcol() const { return m_nResultColumn; }
        void resultcol(int r) { m_nResultColumn=r; }

        unsigned long maxtime() const { return m_MaxTime; }
        void maxtime(unsigned long m) { m_MaxTime=m; }

        double accuracy() const { return m_Accuracy; }
        void accuracy(double a) { m_Accuracy=a; }

        void Run();

    private:

        struct Runner
        {
           Runner(VirtualExperiment *p):pOwner(p) {}
           double operator()(VariablesHolder& v);

           VirtualExperiment *pOwner;
        };

        friend class Runner;

        double getSSRD(std::vector<std::pair<int,double> >& d);
        std::string m_strModelName;
        ObjRef<iface::cellml_api::Model> m_Model;
        int m_nResultColumn;
        typedef std::map<std::wstring,double> PARAMS;
        typedef std::pair<double,double> POINT;
        typedef std::vector<POINT> TIMEPOINTS;

        PARAMS m_Parameters;
        TIMEPOINTS m_Timepoints;
        double m_ReportStep;
        unsigned long m_MaxTime;
        double m_Accuracy;
};


//Virtual experiment group handler
class VEGroup
{
    private:
        VEGroup();
        ~VEGroup();

    public:
        static VEGroup& instance();
        double Evaluate(VariablesHolder& v);
        void add(VirtualExperiment *p);

    protected:
        typedef std::vector<VirtualExperiment *> VE;
        VE experiments;
};

#endif

