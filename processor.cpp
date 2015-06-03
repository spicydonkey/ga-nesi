#include <stdlib.h>
#include <stdio.h>
#include "processor.h"


// The necessary computation based on a workitem's data
	// user must customise this for his needs
double do_compute(std::vector<double>& vals);

// Singleton support to ensure the only instance of the process schedule
Processor& Processor::instance()
{
	static Processor *pInstance=new Processor;

	return *pInstance;
}

// Processor constructor
Processor::Processor()
{
	// Initialisation of Processor object
}

Processor::~Processor()
{
}

// Add a work item to the back of the schedule
void Processor::push(WorkItem* item)
{
	witems.push_back(item);
}

// Remove the workitem in schedule with given key
void Processor::remove_key(int key)
{
	for(WORKITEMS::iterator it=witems.begin();it!=witems.end();)
	{
		if((*it)->key==key)
			it=witems.erase(it);
		else
			++it;
	}
}

// Number of workitems scheduled
int Processor::count()
{
	return witems.size();
}

// Process scheduled workitems
// calls OBESRVER o for each process
// p is a context passed to observer and is transparent to the processor
void Processor::process(OBSERVER o, void *p)
{
	while(witems.size())
	{
		// get next workitem for processing *FIFO*
		WorkItem *pWorkitem=witems.front();
		witems.pop_front();

		// do computation for this workitem
		double answer;
		answer=do_compute(pWorkitem->data);	// compute this workitem's data with user-defined method
		// call observer
		o(pWorkitem,answer,p);
	}
}
