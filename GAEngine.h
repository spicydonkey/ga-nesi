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
		// constructors & destructors
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

		// assignment operator
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

		// allele
        double allele(const std::wstring& name)
        {
            ALLELE::iterator it=find_if(m_Alleles.begin(),m_Alleles.end(),
                   bind1st(pair_equal_to<std::wstring,double>(),name));		// pair_equal_to argument std::wstring& ?
            return (it==m_Alleles.end()?double(0.0):it->second);			// returns 0.0 if name not found and m_Alleles' second value if found 
        }
        double allele(int index)
        {
            return ((index>=0 && index<m_Alleles.size())?m_Alleles[index].second:0.0);	// second val of ith pair of m_Alleles; if index out of range, 0.0
        }
        void allele(const std::wstring& name,double val)		// if m_Alleles has an element with first==name, assign the second to val; push_back <name,val> into m_Allele
        {
            ALLELE::iterator it=find_if(m_Alleles.begin(),m_Alleles.end(),
                   bind1st(pair_equal_to<std::wstring,double>(),name));
            if(it!=m_Alleles.end())
               it->second=val;
            else
               m_Alleles.push_back(std::make_pair<std::wstring,double>(std::wstring(name),double(val)));
        }
        void allele(int index,double val)		// assigns the second memb of the i(ndex)th elem of m_Alleles to val, if in range
        {
            if(index>=0 && index<m_Alleles.size())
               m_Alleles[index].second=val;
        }

		// valid
        bool valid() const { return m_Valid; }
        void valid(bool b) { m_Valid=b; }		// assign m_Valid member to b

		// name
        std::wstring name(int index)
        {
			// rts first memb of ith elem of m_Alleles; if index out of range, an empty string
            return ((index>=0 && index<m_Alleles.size())?m_Alleles[index].first:std::wstring());
        }

		// size
        int size()
        {	
			// rts the size of m_Alleles memb vector
            return m_Alleles.size();
        }

		// [] operator (indexing)
        std::pair<std::wstring,double>& operator[](int index)
        {	
			// pushback <emptystr, 0.0> to m_Alleles until vector index is in range, then return that elem of m_Alleles
            while(m_Alleles.size()<=index)
				m_Alleles.push_back(std::make_pair(std::wstring(),double(0.0)));
			return m_Alleles[index];
        }

		// fitness
        double fitness() const { return m_Fitness; }
        void fitness(double v) { m_Fitness=v; m_Valid=(v!=INFINITY); }		// assigns v to m_Fitness, m_Valid set to finity of v
        
		// </>/== operators (comparator) : Compares the fitness() memb fun of Genomes
		bool operator<(const Genome& other) const
        {
			// returns fitness() < other.fitness(); if at least one invalid, returns this->valid()
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

		// same
        bool same(const Genome& other) const
        {
            /*if(!valid() & other.valid())		// ??? bitwise & on bools
                return false;
            */
			if(!valid() && other.valid())
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

		// var
		// store the alleles data of a genome in a variables holder
			// 'update' all alleles of this genome to v
			// if this genome is missing some alleles, v may contain some alleles unknown to genome
        void var(VariablesHolder& v)
        {
			// iterate through the alleles in this genome
            for(ALLELE::iterator it=m_Alleles.begin();it!=m_Alleles.end();++it)
            {
				// update each allele in v
                v(it->first,it->second);
            }
        }

		// set
        void set(VariablesHolder& v)
        {
			// rebuild alleles from that stored in a varholder

            m_Alleles.clear();	// clear m_Alleles vector
			// m_Alleles is an empty vector of pair

            for(int k=0;;k++)
            {
				std::wstring name=v.name(k);
				if(name.empty())	// k reached the end of m_Var in v
					break;

                m_Alleles.push_back(std::make_pair<std::wstring,double>(name,v(name)));		// append a pair made from m_Vars of v to m_Alleles
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

		// default GA Engine constructor
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

		// Set the maximum population size of GA and resize the population Genome vector accordingly
        void set_borders(int max_population)
        {
            m_MaxPopulation=max_population;
            m_Population.resize(max_population);
        }


		// Initialise the population with randomly generated genomes
        bool Initialise()
        {
			// Create a representative genome for the population
            Genome v;
			
			// exit exceptions
            if(!m_Population.size() || !m_AlleleList.size())
                return false;

			// Attach alleles to the representative genome
			for(typename std::vector<std::wstring>::iterator it=m_AlleleList.begin();it!=m_AlleleList.end();++it)	// iterate through the list of allele names
			{
				// Create each allele with all values initialised to 0.0
				v.allele(*it,(double)0.0);
			}

			// Fill the population with each mutation of the genome
            for(int i=0;i<m_Population.size();i++)
            {
                mutate(std::wstring(),v,true);	// mutate all the alleles of v
                m_Population[i]=v;				// a randomly generated genome enters the population
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
		
		// var_template
		// store in a VarHolder the template allele name list
        void var_template(VariablesHolder& v)
        {
            for(std::vector<std::wstring>::iterator it=m_AlleleList.begin();it!=m_AlleleList.end();++it)
            {
                v(*it,0.0);
            }
        }

		// GetBest
        double GetBest(VariablesHolder& v)
        {
            v=m_bestVariables;
            return m_bestFitness;
        }

		// var_to_workitem
		// returns a pointer to a newly initiated WorkItem that has collated h's allele values
        WorkItem *var_to_workitem(VariablesHolder& h)
        {
            WorkItem *w=new WorkItem;
            w->key=0;
            h.collate(w->data);
            return w;
        }

		// process_workitem
		// assigns the key-th Genome of m_Pop's m_fitness to answer, and deletes the WorkItem
        void process_workitem(WorkItem *w,double answer)
        {
            if(w->key<m_Population.size())
            {
                Genome& g=m_Population[w->key];
                g.fitness(answer);
            }
            delete w;
        }

		// RunGenerations
		/**
		 *	run the GA engine for given number of generations
		 *
		 *
		 **/
        void RunGenerations(int gener)
        {
            VariablesHolder v;
       
            m_Generations=gener;

            //Create initial fitness set
			for(int i=0;i<m_Population.size();i++)
			{
				Genome& g=m_Population[i];	// get the ith Genome in population

				// 'update' all alleles of this genome into temporary variable storage
				g.var(v);	// strange behaviour when this genome is incomplete
				WorkItem *w=var_to_workitem(v);		// initiate a ptr workitem that has collated this genome's allele values
				w->key=i;	// assign this workitem's key accordingly
				Distributor::instance().push(w);	// push this work into the singleton distributor
			}

			// Process the works
			Distributor::instance().process(observer,this);		//?? this must be assigning the fitness of each genome in population: HINT this meaning the GAEngine obj

			std::sort(m_Population.begin(),m_Population.end(),reverse_compare);		// sort m_Population (vector<Genome>) by reverse_compare: in ascending order of fitness
            
			// check if best ftns is assigned. if not, assign it	(fitness minimisation!)
			if(!m_bBestFitnessAssigned || m_bestFitness>m_Population[0].fitness())
            {
                m_bestFitness=m_Population[0].fitness();	// assign min ftns as the best fitness
                m_Population[0].var(m_bestVariables);		// update the bestVars
                m_bBestFitnessAssigned=true;
            }

            print_stage(-1);		// -1 for initial generation

            for(int g=0;g<gener;g++)
            {
				//Do the genetics
				int limit=m_Population.size();
                //Create new population
                POPULATION prev(m_Population);		// prev is a copy of current pop vector
 
                m_Population.clear();		// clear current population vector
                double l=(double)prev.size()-0.5;

				// SELECTION
				// select Genomes from prev gen to carry on; store in m_Population vector
                for(int i=0;i<limit;i++)
                {
                    int mem=select_weighted(prev);		// mem is the randomly selected Genome's index (p.size()-1 when err)
                    //printf("Adding %d to population\n",mem);
                    m_Population.push_back(prev[mem]);	// append the selected Genome into curr Pop vector
                }


                //Do the crossovers
                if(m_crossPartition)
                {
                    std::vector<int> sample;	// initialise an integer vector

                    if(!m_UseBlockSample)
                        build_rnd_sample_rnd(sample,m_CrossProbability*100.0,true);		// fill sample with indices to perform crossover
                    else
                        build_rnd_sample(sample,m_crossPartition,true,true);			// disallow duplicates in building sample (size m_crossPartition)

                    for(int i=0;i<sample.size();i++)
                    {
                        std::vector<int> arena; //initialise arena for breeding
                        
                        arena.push_back(sample[i]);	//ith sample enters arena 
						//bulid tournament sample
                        build_rnd_sample(arena,1,true,true); //another sample enters arena, avoid self for crossbreeding

						//cross the Genomes in arena at a randomly selected crosspoint
	    				cross(m_Population[arena[0]],m_Population[arena[1]],
							(int)rnd_generate(1.0,m_Population[sample[i]].size()));		//crosspoint in [1,allele length of ith sample genome]

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
                        build_rnd_sample_rnd(sample,m_MutationProbability*100.0,false);	//sample vector includes even invalid Genomes
                    else
                        build_rnd_sample(sample,m_mutatePartition,false,false); //allow duplicates and invalid genomes to build sample (size m_mutatePartition)

					//Treatment of invalid genomes in the population
                    for(int i=0;i<m_Population.size();i++)
                    {
						//add all unselected invalid genomes into sample
						if(!m_Population[i].valid() && std::find(sample.begin(),sample.end(),i)==sample.end())
							sample.push_back(i);
                    }

					//Mutate invalid population members
                    for(int i=0;i<sample.size();i++)
                    {
	    				mutate(std::wstring(),m_Population[sample[i]],!(m_Population[sample[i]].valid()));	// mutate-all iff genome is invalid. else mutate approx 1 allele
						m_Population[sample[i]].var(v);
                        Distributor::instance().remove_key(sample[i]); //remove previously requested processing
						WorkItem *w=var_to_workitem(v);
                        m_Population[sample[i]].set(v);		//??
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
				
				// update best fitness
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
			//verbose summary of GA: print all chromosomes of curr gen
            if(verbosity>1)
            {
                printf("--------------------------------------------------------\n");
                for(int j=0;j<m_Population.size();j++)
				{
                    VariablesHolder v;
					
					//print validity, generation #, and fitness of each chromosome
                    printf("%s[%d](%lf) ",(m_Population[j].valid()?" ":"*"),g+1,m_Population[j].fitness());

					//print each chromosome's alleles (name and value)
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
             
			//shorter summary of GA: print currently fittest chromosome
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

		// mutate
        void mutate(const std::wstring& name,Genome& g,bool mutate_all=false)
        {
            double prob=(mutate_all?101.0:100.0/g.size());

            for(int i=0;i<g.size();i++)
            {
                double p=rnd_generate(0.0,100.0);
                
				// mutate if ( p <= prob )
                if(p>prob)		// chance to skip mutation
                   continue;

				/*
				 *	Mutation rate (%):
				 *		100					if mutate_all == true		i.e. mutate all allele
				 *		100.0/g.size()		if mutate_all == false		i.e. mutate approx. just 1 allele		
				 */


				// if		name is non-emtpy wstring, only the allele with matching name may be mutated
				// else if	name is an empty wstring, mutate all alleles
                if(!name.size() || g.name(i)==name)		
                {
                    LIMITS::iterator it=m_Limits.find(g.name(i));	// check for param limits of this allele
                    double val;
                    if(it==m_Limits.end())
                    {
                        //no limits, just use [-RAND_MAX/2,RAND_MAX/2] as a limit
                        val=rnd_generate(-RAND_MAX*0.5,RAND_MAX*0.5);
                    }
                    else
                    {
						// restrict RNG to set limits
                        val=rnd_generate(it->second.first,it->second.second);
                    }
                    g.allele(i,val);	// set the RNG value to allele

                    if(name.size())		// stop mutation if name is non-emtpy (only mutated one allele)
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

		// cross
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

			// check if genomes' alleles have same size and crosspoint lies in valid range
            if(one.size()!=two.size() || one.size()<crosspoint+1)
                return false;

			// genomes same size and valid crosspoint
			// swap the alleles before crosspoint
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

		// build_rnd_sample
        void build_rnd_sample(std::vector<int>& sample,int count,bool reject_duplicates,bool check_valid)
        {
            double limit=(double)m_Population.size()-0.5;

			// append "count" number of randomly selected integers (index for m_Pop) onto sample
            for(;count>0;count--)
            {
                int v;

				// do-while loop: randomly assign an int to v 
					// if reject_duplicates true, build_rnd_sample will not add duplicates to sample
                do
                {
					v=(int)(rnd_generate(0.0,limit));	// v in [0, m_Pop.size()-1]
                    if(check_valid && !m_Population[v].valid())
						continue;		//??if check_valid true, loop until v is a valid Genome? 
						//(BUG - invalid genomes will still be pushed back onto sample)
                }
                while(reject_duplicates && std::find(sample.begin(),sample.end(),v)!=sample.end());
                //Found next value
                sample.push_back(v);	// push onto sample vector
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

		// build_rnd_sample_rnd
        void build_rnd_sample_rnd(std::vector<int>& sample,double prob,bool check_valid)
        {
			// if check_valid, only appends indices of m_Population for which Genomes are valid, at given probability (%)
			// if false (check_valid), appends Genomes at given rate (%)
            for(int i=0;i<m_Population.size();i++)
            {
                if((!check_valid || m_Population[i].valid()) && prob>=rnd_generate(0.0,100.0))
                    sample.push_back(i);
            }
        }

		// select_weighted
        int select_weighted(POPULATION& p)
        {
            double limit=(double)p.size()-0.5;
            double sum=0.0;

			// total sum of population's fitness	(sum is inf if p[i].fitness == 0 or p[i] invalid)
            for(int i=0;i<p.size();i++)
                sum+=(p[i].valid()?1.0/(p[i].fitness()?p[i].fitness():0.000000000001):99999999999.99999);		// !!!shouldn't summand be 0 if p[i] invalid? (BUG?)
				// sum+=(p[i].valid()?1.0/(p[i].fitness()?p[i].fitness():0.000000000001):0.0);

			// use a randomly selected threshold for cum-sum to select i
            double choice=sum*rnd_generate(0.0,1.0);
            for(int i=0;i<p.size();i++)
            {
                choice-=1.0/(p[i].fitness()?p[i].fitness():0.000000000001);
                if(choice<=0.0)
					return i;
            }
			// choice is larger than cum-sum
            return p.size()-1;
        }
};

#endif

