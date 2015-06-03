//Processor class for sequential process scheduling of independent tasks
//Inherent portability for HPC using Distributor class
#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <vector>
#include <list>

//Work item - holds information about
//data to be passed to a compute task
struct WorkItem
{
    int key; //context-dependent value, passed to the observer
    int context; //distribution context
    std::vector<double> data; //data to be distributed
};

/**
 *	Class to schedule independent jobs sequentially as opposed to distributing
 *	Dual class to Distributor
 *	Singleton, only available through instance()
 *	Processor collects workitems to be processed FIFO until process() is called
 *	process then goes through the "ranks" filing workitems to them
 **/
class Processor
{
	private:
		Processor();
		~Processor();
	public:
		typedef bool (*OBSERVER)(WorkItem *, double answer, void *);	// boolean observer function to be called for each returned result

		static Processor& instance();
		void push(WorkItem* item);	// add new workitem for processing
		void remove_key(int key);	// remove all requests with the specified key
		int count();	// number of workitmes
		void process(OBSERVER o, void *ENGINE);	// process workitems calling observer o for each result

	protected:
		typedef std::list<WorkItem*> WORKITEMS;		// list of ptrs to workitem
		WORKITEMS witems;
};


#endif
