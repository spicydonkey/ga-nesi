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


// COMP_FUNC is a function object class for <= comparisons on doubles
#define COMP_FUNC std::less_equal<double>

// define ALLELE as vector of <wstring, double> pairs ( i.e. ALLELE is a vector of 'allele's (pair<wstr,doub>) )
typedef std::vector<std::pair<std::wstring,double> > ALLELE;


/**
 *	VariablesHolder
 *	
 *	+ m_Vars : a storage for ALLELE variable
 *	
 *	
 **/
class VariablesHolder
{
	private:
		ALLELE m_Vars;	//the member that stores alleles in a chromosome form ( i.e. vector <allele=pair<wstring allele_name,double allele_value> > )

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

		//Update an allele (pair<str name, doub value>) in VarHold's m_Var obj and return updated allele value
		double operator()(const std::wstring& name,double val)
		{
			//find if matching allele already exists in VarHold
			ALLELE::iterator it=find_if(m_Vars.begin(),m_Vars.end(),
			   bind1st(pair_equal_to<std::wstring,double>(),name));

			if(it!=m_Vars.end())
			   it->second=val;	//update the allele if it already exists
			else
			   m_Vars.push_back(std::make_pair<std::wstring,double>(std::wstring(name),double(val)));	//add the allele into m_Var if not stored yet
	    
			return val;	//return updated allele value
		}

		//Index search allele name
		std::wstring name(int index)
		{
			//return the name of allele at the index location of allele vector m_Vars (nullwstr if index out of range)
			return ((index>=0 && index<m_Vars.size())?m_Vars[index].first:std::wstring());
		}

		bool exists(const std::wstring& name)
		{
			//return existence of allele of given name in m_Vars vector
			ALLELE::iterator it=find_if(m_Vars.begin(),m_Vars.end(),
			   bind1st(pair_equal_to<std::wstring,double>(),name));
			return (it!=m_Vars.end());
		}

		//Size of a VariablesHolder object
		size_t size()
		{
			//is the size of m_Vars vector i.e. number of alleles in stored in the m_Vars vector
			return m_Vars.size();
		}

		//collate
		//Pushback all allele values held in m_Vars onto a supplied ref to vect<doub>
		void collate(std::vector<double>& v)
		{
			for(ALLELE::iterator it=m_Vars.begin();it!=m_Vars.end();++it)
			{
				v.push_back(it->second);
			}
		}

		//print all alleles held in m_Vars
		void print()
		{
			for(ALLELE::iterator it=m_Vars.begin();it!=m_Vars.end();++it)
			{
				printf("%s->%lf\n",convert(it->first).c_str(),it->second);
			}
		}

		//Fill-up the allele values of m_Vars with a supplied reference to vect<doub>, v
		//returns true iff fillup is executed properly
		bool fillup(std::vector<double>& v)
		{
			//check if the sizes are equal
			if(v.size()!=m_Vars.size())
				return false;

			//assign the second members of each allele in m_Vars
			for(int i=0;i<v.size();i++)
			{
				m_Vars[i].second=v[i];
			}
			return true;
		}
};



/**
 *	VirtualExperiment
 *	
 *	
 *	
 **/
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

		/**
		 *	Runner structure
		 *	
		 *	+ pOwner (ptr to a VirtualExperiment)
		 *
		 **/
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
        
		// Type definitions
		typedef std::map<std::wstring,double>	PARAMS;
        typedef std::pair<double,double>		POINT;
        typedef std::vector<POINT>				TIMEPOINTS;

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
		// get the singleton VE group object
        static VEGroup& instance();

		// TODO
        double Evaluate(VariablesHolder& v);

		// TODO
        void add(VirtualExperiment *p);

    protected:
        typedef std::vector<VirtualExperiment *> VE;
        
		// a vector of pointers to virtual experiments
		VE experiments;
};

#endif

