#ifndef GA_ENGINE_H
#define GA_ENGINE_H
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <limits>
#include <functional>
#include <time.h>
#include <stdlib.h>
#ifdef SUPPORT_MPI
#include <mpi.h>
#endif
#include <limits>
#include "utils.h"
#include "virtexp.h"
#include "distributor.h"
#include <math.h>


extern int verbosity;

//Genome handler
class Genome
{
    private:
        ALLELE m_Alleles;
        double m_Fitness;
        bool m_Valid;

    public:
        Genome():m_Fitness(0.0),m_Valid(true)
        {
        }
        Genome(const Genome& other):m_Fitness(other.m_Fitness),m_Valid(other.m_Valid)
        {
            m_Alleles.assign(other.m_Alleles.begin(),other.m_Alleles.end());
        }
        ~Genome()
        {
        }
        Genome operator=(const Genome& other)
        {
            if(&other!=this)
            {
                m_Alleles.assign(other.m_Alleles.begin(),other.m_Alleles.end());
                m_Fitness=other.m_Fitness;
                m_Valid=other.m_Valid;
             }
            return *this;
        }
        double allele(const std::wstring& name)
        {
            ALLELE::iterator it=find_if(m_Alleles.begin(),m_Alleles.end(),
                   bind1st(pair_equal_to<std::wstring,double>(),name));
            return (it==m_Alleles.end()?double(0.0):it->second);
        }
        double allele(int index)
        {
            return (index>=0 && index<m_Alleles.size()?m_Alleles[index].second:0.0);
        }
        void allele(const std::wstring& name,double val)
        {
            ALLELE::iterator it=find_if(m_Alleles.begin(),m_Alleles.end(),
                   bind1st(pair_equal_to<std::wstring,double>(),name));
            if(it!=m_Alleles.end())
               it->second=val;
            else
               m_Alleles.push_back(std::make_pair<std::wstring,double>(std::wstring(name),double(val)));
        }
        void allele(int index,double val)
        {
            if(index>=0 && index<m_Alleles.size())
               m_Alleles[index].second=val;
        }
        bool valid() const { return m_Valid; }
        void valid(bool b) { m_Valid=b; }

        std::wstring name(int index)
        {
            return ((index>=0 && index<m_Alleles.size())?m_Alleles[index].first:std::wstring());
        }
        int size()
        {
            return m_Alleles.size();
        }
        std::pair<std::wstring,double>& operator[](int index)
        {
             while(m_Alleles.size()<=index)
                m_Alleles.push_back(std::make_pair(std::wstring(),double(0.0)));
             return m_Alleles[index];
        }
        double fitness() const { return m_Fitness; }
        void fitness(double v) { m_Fitness=v; m_Valid=(v!=INFINITY); }
        bool operator<(const Genome& other) const
        {
            if(valid() && other.valid())
                return fitness()<other.fitness();
            return valid();
        }
        bool operator>(const Genome& other) const
        {
            if(valid() && other.valid())
                return fitness()>other.fitness();
            return valid();
        }
        bool operator==(const Genome& other) const
        {
            if(valid() && other.valid())
                return fitness()==other.fitness();
            return false;
        }
        bool same(const Genome& other) const
        {
            if(!valid() & other.valid())
                return false;
            if(other.m_Alleles.size()!=m_Alleles.size())
                return false;
            for(int i=0;i<m_Alleles.size();i++)
            {
                if(m_Alleles[i].first!=other.m_Alleles[i].first ||
                   m_Alleles[i].second!=other.m_Alleles[i].second)
                   return false;
            }
            return true;
        }
        void var(VariablesHolder& v)
        {
            for(ALLELE::iterator it=m_Alleles.begin();it!=m_Alleles.end();++it)
            {
                v(it->first,it->second);
            }
        }
        void set(VariablesHolder& v)
        {
            m_Alleles.clear();

            for(int k=0;;k++)
            {
	        std::wstring name=v.name(k);
	        if(name.empty())
	           break;

                m_Alleles.push_back(std::make_pair<std::wstring,double>(name,v(name)));
            }
        }
};

bool reverse_compare(const Genome& v1,const Genome& v2) { return (v1<v2); }
extern bool observer(WorkItem *w,double answer,void *);

template<typename COMP>
class GAEngine
{
    private:
        typedef std::vector<Genome> POPULATION;
        POPULATION m_Population;
        int m_MaxPopulation;
        std::vector<std::wstring> m_AlleleList;
        double m_CrossProbability;
        double m_MutationProbability;
        int m_crossPartition;
        int m_mutatePartition;
        int m_Generations;
        double m_bestFitness;
        bool m_bBestFitnessAssigned;
        VariablesHolder m_bestVariables;
        bool m_UseBlockSample;

    public:
        typedef Genome GENOME;

        GAEngine():m_MaxPopulation(0),m_Generations(1),
                   m_CrossProbability(0.2),m_MutationProbability(0.01),
                   m_bBestFitnessAssigned(false),m_UseBlockSample(false),
                   m_crossPartition(0),m_mutatePartition(0)
        {
        }
        ~GAEngine()
        {
        }

        double& prob_cross() { return m_CrossProbability; }
        double& prob_mutate() { return m_MutationProbability; }
        int& part_cross() { return m_crossPartition; }
        int& part_mutate() { return m_mutatePartition; }

        bool& block_sample() { return m_UseBlockSample; }

        void set_borders(int max_population)
        {
            m_MaxPopulation=max_population;
            m_Population.resize(max_population);
        }

        bool Initialise()
        {
            Genome v;

            
            if(!m_Population.size() || !m_AlleleList.size())
                return false;

	    for(typename std::vector<std::wstring>::iterator it=m_AlleleList.begin();it!=m_AlleleList.end();++it)
	    {
                v.allele(*it,(double)0.0);
	    }
            for(int i=0;i<m_Population.size();i++)
            {
                mutate(std::wstring(),v,true);
                m_Population[i]=v;
            }
            return true;
        }
        int size() { return m_Population.size(); }
        void AddAllele(const std::wstring& name)
        {
            m_AlleleList.push_back(name);
        }
        void AddLimit(const std::wstring& name,double lower,double upper)
        {
            m_Limits[name]=std::make_pair(double(lower),double(upper));
        }
        void var_template(VariablesHolder& v)
        {
            for(std::vector<std::wstring>::iterator it=m_AlleleList.begin();it!=m_AlleleList.end();++it)
            {
                v(*it,0.0);
            }
        }

        double GetBest(VariablesHolder& v)
        {
            v=m_bestVariables;
            return m_bestFitness;
        }

        WorkItem *var_to_workitem(VariablesHolder& h)
        {
            WorkItem *w=new WorkItem;
            w->key=0;
            h.collate(w->data);
            return w;
        }

        void process_workitem(WorkItem *w,double answer)
        {
            if(w->key<m_Population.size())
            {
                Genome& g=m_Population[w->key];
                g.fitness(answer);
            }
            delete w;
        }

        void RunGenerations(int gener)
        {
            VariablesHolder v;
       
            m_Generations=gener;

            //Create initial fitness set
	    for(int i=0;i<m_Population.size();i++)
	    {
		Genome& g=m_Population[i];

		g.var(v);
		WorkItem *w=var_to_workitem(v);
		w->key=i;
		Distributor::instance().push(w);
	    }
	    Distributor::instance().process(observer,this);
	    std::sort(m_Population.begin(),m_Population.end(),reverse_compare);
            if(!m_bBestFitnessAssigned || m_bestFitness>m_Population[0].fitness())
            {
                m_bestFitness=m_Population[0].fitness();
                m_Population[0].var(m_bestVariables);
                m_bBestFitnessAssigned=true;
            }

            print_stage(-1);

            for(int g=0;g<gener;g++)
            {
		//Do the genetics
	        int limit=m_Population.size();
                //Create new population
                POPULATION prev(m_Population);
 
                m_Population.clear();
                double l=(double)prev.size()-0.5;
                for(int i=0;i<limit;i++)
                {
                    int mem=select_weighted(prev);
                    //printf("Adding %d to population\n",mem);
                    m_Population.push_back(prev[mem]);
                }


                //Do the crossovers
                if(m_crossPartition)
                {
                    std::vector<int> sample;

                    if(!m_UseBlockSample)
                        build_rnd_sample_rnd(sample,m_CrossProbability*100.0,true);
                    else
                        build_rnd_sample(sample,m_crossPartition,true,true); //disallow duplicates
                    for(int i=0;i<sample.size();i++)
                    {
                        std::vector<int> arena; //arena for breeding
                        
                        arena.push_back(sample[i]); //avoid self for crossbreeding
                        build_rnd_sample(arena,1,true,true); //build tournament sample
	    	        cross(m_Population[arena[0]],m_Population[arena[1]],
                              (int)rnd_generate(1.0,m_Population[sample[i]].size()));

                        for(int j=0;j<2;j++)
                        {
			    m_Population[arena[j]].var(v);
                            Distributor::instance().remove_key(arena[j]); //remove previously requested processing
			    WorkItem *w=var_to_workitem(v);
			    w->key=arena[j];
			    Distributor::instance().push(w);
                        }
                    }
                }
                //Do the mutations
                if(m_mutatePartition)
                {
                    std::vector<int> sample;

                    if(!m_UseBlockSample)
                        build_rnd_sample_rnd(sample,m_MutationProbability*100.0,false);
                    else
                        build_rnd_sample(sample,m_mutatePartition,false,false); //allow duplicates

                    //Mutate invalid population members
                    for(int i=0;i<m_Population.size();i++)
                    {
                         if(!m_Population[i].valid() && std::find(sample.begin(),sample.end(),i)==sample.end())
                             sample.push_back(i);
                    }
                    for(int i=0;i<sample.size();i++)
                    {
	    	        mutate(std::wstring(),m_Population[sample[i]],!m_Population[sample[i]].valid());
			m_Population[sample[i]].var(v);
                        Distributor::instance().remove_key(sample[i]); //remove previously requested processing
			WorkItem *w=var_to_workitem(v);
                        m_Population[sample[i]].set(v);
			w->key=sample[i];
			Distributor::instance().push(w);
                    }
                }
		//Run the distribution
		Distributor::instance().process(observer,this);
		std::sort(m_Population.begin(),m_Population.end(),reverse_compare);

                if(m_Population.size()>m_MaxPopulation)
                {
                      //Cull it
                      m_Population.erase(m_Population.begin()+m_MaxPopulation,m_Population.end());
                }
 
                if(!m_bBestFitnessAssigned || m_bestFitness>m_Population[0].fitness())
                {
                    m_bestFitness=m_Population[0].fitness();
                    m_Population[0].var(m_bestVariables);
                    m_bBestFitnessAssigned=true;
                }
                print_stage(g);
            }
        }

    private:
        typedef std::map<std::wstring,std::pair<double,double> > LIMITS;
        LIMITS m_Limits;

        void print_stage(int g)
        {
            if(verbosity>1)
            {
                printf("--------------------------------------------------------\n");
                for(int j=0;j<m_Population.size();j++)
	        {
                    VariablesHolder v;
 
                    printf("%s[%d](%lf) ",(m_Population[j].valid()?" ":"*"),g+1,m_Population[j].fitness());
                    for(int k=0;;k++)
                    {
                        m_Population[j].var(v);
	                std::wstring name=v.name(k);
                               
	    	        if(name.empty())
			    break;
			printf("%s=%lf   ",convert(name).c_str(),v(name));
                    }
                    printf("\n");
                }
                printf("--------------------------------------------------------\n");
             } 
             else if(verbosity==1)
             {
                 VariablesHolder v;

                 double f=GetBest(v);
	         printf("Generation %d. Best fitness: %lf\n",g+1,f);
		 for(int k=0;;k++)
		 {
		     std::wstring name=v.name(k);
		     if(name.empty())
		         break;
		     printf("%s=%lf    ",convert(name).c_str(),v(name));
		 }
                 printf("\n");
                 printf("--------------------------------------------------------\n");
              }
        }

        void mutate(const std::wstring& name,Genome& g,bool mutate_all=false)
        {
            double prob=(mutate_all?101.0:100.0/g.size());

            for(int i=0;i<g.size();i++)
            {
                double p=rnd_generate(0.0,100.0);
                
                if(p>prob)
                   continue;
                if(!name.size() || g.name(i)==name)
                {
                    LIMITS::iterator it=m_Limits.find(g.name(i));
                    double val;
                    if(it==m_Limits.end())
                    {
                        //no limits, just use [-RAND_MAX/2,RAND_MAX/2] as a limit
                        val=rnd_generate(-RAND_MAX*0.5,RAND_MAX*0.5);
                    }
                    else
                    {
                        val=rnd_generate(it->second.first,it->second.second);
                    }
                    g.allele(i,val);

                    if(name.size())
                        break;
                }
            }
        }
        void mutate(double probability,Genome& g,int count=-1)
        {
            int cnt=0;
            double prob=rnd_generate(0.0,100.0);
 
            for(int i=0;i<g.size();i++)
            {
                if(prob<=probability)
                {
                    LIMITS::iterator it=m_Limits.find(g.name(i));
                    double val;
                    if(it==m_Limits.end())
                    {
                        //no limits, just use [-RAND_MAX/2,RAND_MAX] as a limit
                        val=rnd_generate(-RAND_MAX*0.5,RAND_MAX*0.5);
                    }
                    else
                    {
                        val=rnd_generate(it->second.first,it->second.second);
                    }
                    g.allele(i,val);
                    cnt++;
                    if(count>=0 && cnt>=count)
                       break;
                }
            }
        }
        bool cross(Genome& one,Genome& two,int crosspoint,Genome& out)
        {
            if(one.size()!=two.size() || one.size()<crosspoint+1)
                return false;
            for(int i=0;i<crosspoint;i++)
                out[i]=one[i];
            for(int i=crosspoint;i<one.size();i++)
                out[i]=two[i];
            return true;
        }
        bool cross(Genome& one,Genome& two,int crosspoint)
        {
            Genome n1,n2;

            if(one.size()!=two.size() || one.size()<crosspoint+1)
                return false;

            for(int i=0;i<crosspoint;i++)
            {
                n1[i]=two[i];
                n2[i]=one[i];
            }
            
            for(int i=crosspoint;i<one.size();i++)
            {
                n2[i]=two[i];
                n1[i]=one[i];
            }

            one=n1;
            two=n2;

            return true;
        }
        void build_rnd_sample(std::vector<int>& sample,int count,bool reject_duplicates,bool check_valid)
        {
            double limit=(double)m_Population.size()-0.5;

            for(;count>0;count--)
            {
                int v;

                do
                {
                    v=(int)(rnd_generate(0.0,limit));
                    if(check_valid && !m_Population[v].valid())
                       continue;
                }
                while(reject_duplicates && std::find(sample.begin(),sample.end(),v)!=sample.end());
                //Found next value
                sample.push_back(v);
            }
        }
        void build_rnd_sample_tournament(std::vector<int>& sample,int count,bool reject_duplicates,bool check_valid)
        {
            double limit=(double)m_Population.size()-0.5;
            count*=2; //create tournament pairs
            int index=sample.size();

            for(;count>0;count--)
            {
                int v;

                do
                {
                    v=(int)(rnd_generate(0.0,limit));
                    if(check_valid && !m_Population[v].valid())
                       continue;
                }
                while(reject_duplicates && std::find(sample.begin(),sample.end(),v)!=sample.end());
                //Found next value
                sample.push_back(v);
            }
            //let the fight begins!
            for(int i=index;i<sample.size();i++)
            {
                Genome& one=m_Population[sample[i]];    
                Genome& two=m_Population[sample[i+1]];
                
                if(one>two)
                    sample.erase(sample.begin()+i);
                else  
                    sample.erase(sample.begin()+i+1);
            }     
        }
        void build_rnd_sample_rnd(std::vector<int>& sample,double prob,bool check_valid)
        {

            for(int i=0;i<m_Population.size();i++)
            {
                if((!check_valid || m_Population[i].valid()) && prob>=rnd_generate(0.0,100.0))
                    sample.push_back(i);
            }
        }
        int select_weighted(POPULATION& p)
        {
            double limit=(double)p.size()-0.5;
            double sum=0.0;

            for(int i=0;i<p.size();i++)
                sum+=(p[i].valid()?1.0/(p[i].fitness()?p[i].fitness():0.000000000001):99999999999.99999);

            double choice=sum*rnd_generate(0.0,1.0);
            for(int i=0;i<p.size();i++)
            {
                choice-=1.0/(p[i].fitness()?p[i].fitness():0.000000000001);
                if(choice<=0.0)
                     return i;
            }
            return p.size()-1;
        }
};

#endif

