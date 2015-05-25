//Distributor class provides work sharing support
//using MPI for communications
#ifndef DISTRIBUTOR_H
#define DISTRIBUTOR_H

#include <vector>
#include <list>

#define TAG_QUIT 0x100


//Work item - holds information about
//data to be passed to a compute task
struct WorkItem
{
    int key; //context-dependent value, passed to the observer
    int context; //distribution context
    std::vector<double> data; //data to be distributed
};

//Class to handle job distribution via MPI
//Singleton, only available through instance()
//Distributor collects work items to be processed until process() is called
//process then goes through the ranks filing workitems to them
//and calls the OBSERVER callback for every result received
class Distributor
{
    private:
        Distributor();
        ~Distributor();
    public:
        typedef bool (*OBSERVER)(WorkItem *,double answer,void *); //observer function to be called for each returned result

        static Distributor& instance();
        void push(WorkItem* item); //Add new workitem for processing
        void remove_key(int key); //remove all requests with the specified key
        int count(); //number of workitems
        void process(OBSERVER o,void *d); //process workitems calling observer o for each result
        void finish(); //terminate MPI chain, must be called before MPI_Finalize

    protected:
        typedef std::list<WorkItem*> WORKITEMS;
        typedef std::vector<std::pair<bool,WorkItem*> > RANKS;
        WORKITEMS witems;
        RANKS ranks;        
};


#endif

