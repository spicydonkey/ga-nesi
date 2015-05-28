#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "distributor.h"

using namespace std;


double do_compute(std::vector<double>& vals);


//Singleton support to ensure the only instance of the distributor to exist
Distributor& Distributor::instance()
{
    static Distributor *pInstance=new Distributor;
    
    return *pInstance;
}

//Distributor constructor
Distributor::Distributor()
{
    int nproc;
    
    //Initialise ranks array
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    ranks.resize(nproc);
    for(int i=0;i<nproc;i++)
    {
        ranks[i].first=false; //mark all tasks to be not busy
    }
}

Distributor::~Distributor()
{
}


//Add a work item for further processing
void Distributor::push(WorkItem* item)
{
    witems.push_back(item);
}

//returns number of workitems registered
int Distributor::count()
{
    return witems.size();
}

void Distributor::remove_key(int key)
{
    for(WORKITEMS::iterator it=witems.begin();it!=witems.end();)
    {
        if((*it)->key==key)
            it=witems.erase(it);
        else
            ++it;
    }
}

//Process registered workitems
//calls OBSERVER o for each new reply
//p is a context passed to observer and is transparent for the distributor
void Distributor::process(Distributor::OBSERVER o,void *p)
{
    int in_process=0;


    while(witems.size())
    {
        int i=1;

        //get next workitem for processing
        WorkItem *workitem=witems.front();
        witems.pop_front();

        //find a rank to send the workitem to
        for(;i<ranks.size();i++)
            if(!ranks[i].first)
               break; //found available rank
        
        //check if an available rank is found
        if(i<ranks.size())
        {
            ranks[i].first=true; //Found available rank
            ranks[i].second=workitem;
            workitem->context=time(NULL); //save time for adding load balancing later
            //Request processing
            MPI_Send(&workitem->data[0],workitem->data.size(),MPI_DOUBLE,i,0,MPI_COMM_WORLD);
            in_process++;
        }
        else
        {
            //we are the only one available - do compute
            double answer;

            answer=do_compute(workitem->data); // compute this workitem's data, returns the residual
            o(workitem,answer,p); //call observer
            //get data back
            if(in_process)
            {
                MPI_Status stat;
                int r;
                int flag=0;

                //some ranks are still processing, check if anything is returned
                while(true)
                {
                    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&stat);
                    if(flag)
                    {
                        //There is data available
		        r=stat.MPI_SOURCE;
		        MPI_Recv(&answer,1,MPI_DOUBLE,stat.MPI_SOURCE,0,MPI_COMM_WORLD,&stat);
		        o(ranks[r].second,answer,p);                 
		        ranks[r].first=false;
		        in_process--;
                    }
                    else
                       break;
                }
            }
        }
    }
    //all done - wait for the rest if anything left
    while(in_process)
    {
	MPI_Status stat;
	int r;
        double answer;

        //Wait until someting is returned - can be ommitted
	MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	r=stat.MPI_SOURCE;
	MPI_Recv(&answer,1,MPI_DOUBLE,stat.MPI_SOURCE,0,MPI_COMM_WORLD,&stat);
	o(ranks[r].second,answer,p);                 
	ranks[r].first=false;
	in_process--;
    }
}


//finalize processing and notifies all the ranks about
//requested end of service
void Distributor::finish()
{
    //Indicate end-of-work to all the slaves
    for(int i=1;i<ranks.size();i++)
    {
        MPI_Send(&i,1,MPI_INT,i,TAG_QUIT,MPI_COMM_WORLD);
    }
}

