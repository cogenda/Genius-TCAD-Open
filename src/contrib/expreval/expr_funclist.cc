// File:    funclist.cpp
// Author:  Brian Vanderburg II
// Purpose: Function list
//------------------------------------------------------------------------------


// Includes
#include <new>

#include "expr_funclist.h"
#include "expr_except.h"
#include "expr_node.h"

using namespace std;
using namespace ExprEval;

// Function factory
//------------------------------------------------------------------------------

// Constructor
FunctionFactory::FunctionFactory()
    {
    }
    
// Destructor
FunctionFactory::~FunctionFactory()
    {
    }
    
// Create
FunctionNode *FunctionFactory::Create(Expression *expr)
    {
    FunctionNode *n = DoCreate(expr);
    if(n)
        n->m_factory = this;
        
    return n;
    }
        

// Function list
//------------------------------------------------------------------------------

// Constructor
FunctionList::FunctionList()
    {
    }
                    
// Destructor
FunctionList::~FunctionList()
    {
    // Free function factories
    Clear();
    }
            
// Add factory to list
void FunctionList::Add(FunctionFactory *factory)
    {
    // Check it
    if(factory == 0)
        throw(NullPointerException("FunctionList::Add"));
    
    // Make sure it does not exist
    size_type pos;
    
    for(pos  = 0; pos < m_functions.size(); pos++)
        {
        if(m_functions[pos]->GetName() == factory->GetName())
            throw(AlreadyExistsException(factory->GetName()));
        }
    
    m_functions.push_back(factory);
    }
    
// Create a node for a function
FunctionNode *FunctionList::Create(const string &name, Expression *expr)
    {
    // Make sure pointer exists
    if(expr == 0)
        throw(NullPointerException("FunctionList::Create"));
    
    size_type pos;
    
    for(pos = 0; pos < m_functions.size(); pos++)
        {
        if(m_functions[pos]->GetName() == name)
            {
            // Found it
            return m_functions[pos]->Create(expr);
            }
        }
        
    // Not found
    return 0;
    }
            
// FunctionList::AddDefaultFunctions is located in func.cpp
// along with the default function factories            
          
// Free function list
void FunctionList::Clear()
    {
    size_type pos;
    
    for(pos = 0; pos < m_functions.size(); pos++)
        {
        delete m_functions[pos];
        }
    }




