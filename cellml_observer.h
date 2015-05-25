//Class LocalProgressObserver implements
//observer role for CellML API purposes
#ifndef CELLML_OBSERVER_H
#define CELLML_OBSERVER_H
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "cellml-api-cxx-support.hpp"
#include "IfaceCellML_APISPEC.hxx"
#include "CellMLBootstrap.hpp"
#include <string>
#include "utils.h"


class LocalProgressObserver:public iface::cellml_services::IntegrationProgressObserver
{
public:
  LocalProgressObserver(iface::cellml_services::CellMLCompiledModel* aCCM)
    : mRefcount(1), bFinished(false),bFailed(false)
  {
    mCCM = aCCM;
    mCCM->add_ref();
    mCI = mCCM->codeInformation();

    iface::cellml_services::ComputationTargetIterator* cti =
      mCI->iterateTargets();
    while (true)
    {
      iface::cellml_services::ComputationTarget* ct = cti->nextComputationTarget();
      if (ct == NULL)
        break;
      if ((ct->type() == iface::cellml_services::STATE_VARIABLE ||
           ct->type() == iface::cellml_services::ALGEBRAIC ||
           ct->type() == iface::cellml_services::VARIABLE_OF_INTEGRATION) &&
          ct->degree() == 0)
      {
        iface::cellml_api::CellMLVariable* source = ct->variable();
        std::string n = convert(source->name());
        source->release_ref();
      }
      ct->release_ref();
    }
    cti->release_ref();
  }

  ~LocalProgressObserver()
  {
    mCCM->release_ref();
    mCI->release_ref();
  }

  void add_ref()
    throw(std::exception&)
  {
    mRefcount++;
  }

  void release_ref()
    throw(std::exception&)
  {
    mRefcount--;
    if (mRefcount == 0)
      delete this;
  }

  std::string objid()
    throw (std::exception&)
  {
    return std::string("singletonLocalProgressObserver");
  }

  void* query_interface(const std::string& iface)
    throw (std::exception&)
  {
    if (iface=="xpcom::IObject")
      return static_cast< ::iface::XPCOM::IObject* >(this);
    else if (iface=="cellml_services::IntegrationProgressObserver")
      return
        static_cast< ::iface::cellml_services::IntegrationProgressObserver*>
        (this);
    return NULL;
  }

virtual std::vector<std::string> supported_interfaces() throw (std::exception&)
{
    std::vector<std::string> v;

    v.push_back("xpcom::IObject");
    v.push_back("cellml_services::IntegrationProgressObserver");
    return v;
}


  void computedConstants(const std::vector<double>& results)
    throw (std::exception&)
  {
#if 0
    iface::cellml_services::ComputationTargetIterator* cti =
      mCI->iterateTargets();
    while (true)
    {
      iface::cellml_services::ComputationTarget* ct = cti->nextComputationTarget();
      if (ct == NULL)
        break;
      if (ct->type() == iface::cellml_services::CONSTANT &&
          ct->degree() == 0)
      {
        iface::cellml_api::CellMLVariable* source = ct->variable();
        std::string n = convert(source->name());
        source->release_ref();
        printf("# Computed constant: %S = %e\n", n.c_str(), values[ct->assignedIndex()]);
      }
      ct->release_ref();
    }
    cti->release_ref();
#endif
  }

  void results(const std::vector<double>& results)
    throw (std::exception&)
  {
      m_Results.insert(m_Results.end(),results.begin(),results.end());
  }

//Public interface to the observer data
//
//results of the computations are returned in res vector
  int GetResults(std::vector<double>& res)
  {
    uint32_t aic = mCI->algebraicIndexCount();
    uint32_t ric = mCI->rateIndexCount();
    uint32_t recsize = 2 * ric + aic + 1;

    if (recsize == 0)
      return 0;
    res.assign(m_Results.begin(),m_Results.end());
    return recsize;
  }

//Marks computation process as finished
  void done()
    throw (std::exception&)
  {
    bFinished = true;
  }

//Marks computation process as finished
//and prints out error message
  void failed(const std::string& errmsg)
    throw (std::exception&)
  {
    fprintf(stderr,"# Integration failed (%s)\n", errmsg.c_str());
    bFinished = true;
    bFailed=true;
  }

//true if compute is done
  bool finished() const { return bFinished; }
  bool failed() const { return bFailed; }

private:
  iface::cellml_services::CellMLCompiledModel* mCCM;
  iface::cellml_services::CodeInformation* mCI;
  uint32_t mRefcount;
  bool bFinished;
  bool bFailed;
  std::vector<double> m_Results;
};

#endif

