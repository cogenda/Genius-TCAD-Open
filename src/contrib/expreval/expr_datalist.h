// File:    datalist.h
// Author:  Brian Vanderburg II
// Purpose: Data list used for arbitrary data support
//------------------------------------------------------------------------------


#ifndef __EXPREVAL_DATALIST_H
#define __EXPREVAL_DATALIST_H

// Includes
#include <string>
#include <vector>

// Part of expreval namespace
namespace ExprEval
    {
    // Base data type
    //--------------------------------------------------------------------------
    class Data
        {
        public:
            virtual ~Data();
            
            virtual long GetType() const = 0;
            
            static long CreateType(char a, char b, char c, char d);
        };
        
    // Data entry
    //--------------------------------------------------------------------------    
    class DataEntry
        {
        public:
            DataEntry();
            ~DataEntry();
            
            Data *SetData(Data *data);
            Data *GetData();
        
        private:
            Data *m_data;
        };
    
    // Data list item
    //--------------------------------------------------------------------------
    class DataListItem
        {
        public:
            DataListItem(const ::std::string &name, Data *data = NULL);
            
            const ::std::string &GetName() const;
            DataEntry *GetDataEntry();            
            
        private:
            ::std::string m_name; // Name of value
            DataEntry m_data;
        };
        
    // Data list
    //--------------------------------------------------------------------------
    class DataList
        {
        public:
            typedef ::std::vector<DataListItem*>::size_type size_type;
                    
            DataList();
            ~DataList();
            
            // Add data list item
            void Add(const ::std::string &name, Data *data = NULL);
            
            // Get the address of a data entry
            DataEntry *GetDataEntry(const ::std::string &name);
            
            // Enumerate data
            size_type Count() const;
            DataEntry *Item(size_type pos);
            
            // Reset items by deleting data
            void Reset();
            
            // Free items and clear memory
            void Clear();
        
        private:
            ::std::vector<DataListItem*> m_data;
        };
    }

#endif // __EXPREVAL_DATALIST_H

