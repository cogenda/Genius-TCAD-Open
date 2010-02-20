// File:    expr.cc
// Author:  Brian Vanderburg II
// Purpose: Expression object
//------------------------------------------------------------------------------

// Includes
#include <new>
#include <memory>

#include "expr.h"
#include "expr_parser.h"
#include "expr_node.h"
#include "expr_except.h"

using namespace std;
using namespace ExprEval;


// Expression object
//------------------------------------------------------------------------------

// Constructor
Expression::Expression() : m_vlist(0), m_flist(0), m_dlist(0), m_expr(0)
    {
    m_abortcount = 200000;
    m_abortreset = 200000;
    }
    
// Destructor
Expression::~Expression()
    {
    // Delete expression nodes
    delete m_expr;
    }

// Set value list
void Expression::SetValueList(ValueList *vlist)
    {
    m_vlist = vlist;
    }
    
// Get value list
ValueList *Expression::GetValueList() const
    {
    return m_vlist;
    }

// Set function list
void Expression::SetFunctionList(FunctionList *flist)
    {
    m_flist = flist;
    }
    
// Get function list
FunctionList *Expression::GetFunctionList() const
    {
    return m_flist;
    }     
    
// Set data list
void Expression::SetDataList(DataList *dlist)
    {
    m_dlist = dlist;
    }
    
// Get data list
DataList *Expression::GetDataList() const
    {
    return m_dlist;
    }        
            
// Test for an abort
bool Expression::DoTestAbort()
    {
    // Derive a class to test abort
    return false;
    }
    
// Test for an abort
void Expression::TestAbort(bool force)
    {
    if(force)
        {
        // Test for an abort now
        if(DoTestAbort())
            {
            throw(AbortException());
            }
        }
    else
        {
        // Test only if abort count is 0
        if(m_abortcount == 0)
            {
            // Reset count
            m_abortcount = m_abortreset;
            
            // Test abort
            if(DoTestAbort())
                {
                throw(AbortException());
                }
            }
        else
            {
            // Decrease abort count
            m_abortcount--;
            }
        }
    }

// Set test abort count
void Expression::SetTestAbortCount(unsigned long count)
    {
    m_abortreset = count;
    if(m_abortcount > count)
        m_abortcount = count;
    }
            
// Parse expression
void Expression::Parse(const string &exstr)
    {
    // Clear the expression if needed
    if(m_expr)
        Clear();
        
    // Create parser
    auto_ptr<Parser> p(new Parser(this));
    
    // Parse the expression
    m_expr = p->Parse(exstr);
    }
    
// Clear the expression
void Expression::Clear()
    {
    delete m_expr;
    m_expr = 0; 
    }

// Evaluate an expression
double Expression::Evaluate()
    {
    if(m_expr)
        {
        return m_expr->Evaluate();
        }
    else
        {
        throw(EmptyExpressionException());
        }    
    }
            
