// File:    vallist.cpp
// Author:  Brian Vanderburg II
// Purpose: Value list used for variable and constants
//------------------------------------------------------------------------------


// Includes
#include <new>
#include <memory>

#include "expr_defs.h"
#include "expr_vallist.h"
#include "expr_except.h"

using namespace std;
using namespace ExprEval;

// ValueListItem
//------------------------------------------------------------------------------

// Constructor for internal value
ValueListItem::ValueListItem(const string &name, double def, bool constant)
    {
    m_name = name;
    m_constant = constant;
    m_value = m_def = def;
    m_ptr = 0;
    }

// Constructor for external value    
ValueListItem::ValueListItem(const string &name, double *ptr, double def, bool constant)
    {
    m_name = name;
    m_constant = constant;
    m_value = m_def = def;    
    m_ptr = ptr;
    
    if(m_ptr)
        *m_ptr = def;
    else
        throw(NullPointerException("ValueListItem::ValueListItem"));        
    }    
    
// Get the name
const string &ValueListItem::GetName() const
    {
    return m_name;
    }
        
// Return if it is constant
bool ValueListItem::IsConstant() const
    {
    return m_constant;
    }
    
// Get value address
double *ValueListItem::GetAddress()
    {
    return m_ptr ? m_ptr : &m_value;
    } 
    
// Reset to default value
void ValueListItem::Reset()
    {
    if(m_ptr)
        *m_ptr = m_def;
    else
        m_value = m_def;
    }    

    
// ValueList
//------------------------------------------------------------------------------

// Constructor
ValueList::ValueList()
    {
    }
    
// Destructor
ValueList::~ValueList()
    {
    Clear();
    }    
    
// Add value to list
void ValueList::Add(const string &name, double def, bool constant)
    {
    // Ensure value does not already exist
    if(GetAddress(name))
        {
        throw(AlreadyExistsException(name));
        }
    else
        {
        // Create value
        auto_ptr<ValueListItem> i(new ValueListItem(name, def ,constant));
        
        // Add value to list
        m_values.push_back(i.get());
        i.release();
        }
    }
    
// Add an external value to the list
void ValueList::AddAddress(const string &name, double *ptr, double def, bool constant)
    {
    if(GetAddress(name))
        {
        throw(AlreadyExistsException(name));
        }
    else if(ptr == 0)
        {
        throw(NullPointerException("ValueList::AddAddress"));
        }
    else
        {
        // Create value
        auto_ptr<ValueListItem> i(new ValueListItem(name, ptr, def, constant));
        
        // Add value to list
        m_values.push_back(i.get());
        i.release();
        }
    }    
    
// Get the address of the value, internal or external
double *ValueList::GetAddress(const string &name) const
    {
    size_type pos;
    
    for(pos = 0; pos < m_values.size(); pos++)
        {
        if(m_values[pos]->GetName() == name)
            {
            // The name matches
            return m_values[pos]->GetAddress();
            }
        }
        
    // No item found
    return 0;
    }    
    
// Is the value a constant
bool ValueList::IsConstant(const string &name) const
    {
    size_type pos;
    
    for(pos = 0; pos < m_values.size(); pos++)
        {
        if(m_values[pos]->GetName() == name && m_values[pos]->IsConstant())
            {
            return true;
            }
        }
        
    return false;
    }   
    
// Number of values in the list
ValueList::size_type ValueList::Count() const
    {
    return m_values.size();
    }
    
// Get an item
void ValueList::Item(size_type pos, string *name, double *value) const
    {
    if(name)
        *name = m_values[pos]->GetName();
        
    if(value)
        *value = *(m_values[pos]->GetAddress());
    }
    
// Add some default values
void ValueList::AddDefaultValues()
    {
    // Math constant 'e'
    Add("E", EXPREVAL_E, true);
    
    // Math constant PI
    Add("PI", EXPREVAL_PI, true);
    }
    
// Reset values
void ValueList::Reset()
    {
    size_type pos;
    
    for(pos = 0; pos < m_values.size(); pos++)
        {
        m_values[pos]->Reset();
        }
    }
    
// Free values
void ValueList::Clear()
    {
    size_type pos;
    
    for(pos = 0; pos < m_values.size(); pos++)
        {
        delete m_values[pos];
        }
    }        
                                                  
