// File:    datalist.cc
// Author:  Brian Vanderburg II
// Purpose: Data list for arbitrary data
//------------------------------------------------------------------------------


// Includes
#include <new>
#include <memory>

#include "expr_defs.h"
#include "expr_datalist.h"
#include "expr_except.h"

using namespace std;
using namespace ExprEval;

// Data
//------------------------------------------------------------------------------
Data::~Data()
    {
    }
    
long Data::CreateType(char a, char b, char c, char d)
    {
    return (long(a) << 24) | (long(b) << 16) | (long(c) << 8) | long(d);
    }    
    
// DataEntry
//------------------------------------------------------------------------------
DataEntry::DataEntry() : m_data(0)
    {
    }
    
DataEntry::~DataEntry()
    {
    delete m_data;
    }    
    
Data *DataEntry::SetData(Data *data)
    {
    Data *tmp;
    
    tmp = m_data;
    m_data = data;
    
    return tmp;
    }   
    
Data *DataEntry::GetData()
    {
    return m_data;
    }    

// DataListItem
//------------------------------------------------------------------------------

DataListItem::DataListItem(const string &name, Data *data)
    {
    m_name = name;
    m_data.SetData(data);
    }

// Get the name
const string &DataListItem::GetName() const
    {
    return m_name;
    }
        
// Get data entry
DataEntry *DataListItem::GetDataEntry()
    {
    return &m_data;
    } 
    
// DataList
//------------------------------------------------------------------------------

// Constructor
DataList::DataList()
    {
    }
    
// Destructor
DataList::~DataList()
    {
    Clear();
    }    
    
// Add data to list
void DataList::Add(const string &name, Data *data)
    {
    // Ensure data name does not already exist
    if(GetDataEntry(name))
        {
        throw(AlreadyExistsException(name));
        }
    else
        {
        // Create data item
        auto_ptr<DataListItem> i(new DataListItem(name, data));
        
        // Add to list
        m_data.push_back(i.get());
        i.release();
        }
    }
    
// Get the address of the data entry
DataEntry *DataList::GetDataEntry(const string &name)
    {
    size_type pos;
    
    for(pos = 0; pos < m_data.size(); pos++)
        {
        if(m_data[pos]->GetName() == name)
            {
            // The name matches
            return m_data[pos]->GetDataEntry();
            }
        }
        
    // No item found
    return 0;
    }    
    
// Number of data entries
DataList::size_type DataList::Count() const
    {
    return m_data.size();
    }
    
// Get an entry
DataEntry *DataList::Item(size_type pos)
    {
    return m_data[pos]->GetDataEntry();
    }
    
// Reset data
void DataList::Reset()
    {
    size_type pos;
    
    for(pos = 0; pos < m_data.size(); pos++)
        {
        DataEntry *e = m_data[pos]->GetDataEntry();
        
        if(e)
            delete e->SetData(0);
        }
    }
    
// Free data
void DataList::Clear()
    {
    size_type pos;
    
    for(pos = 0; pos < m_data.size(); pos++)
        {
        delete m_data[pos];
        }
    }        
                                                  
