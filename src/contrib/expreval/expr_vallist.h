// File:    vallist.h
// Author:  Brian Vanderburg II
// Purpose: Value list used for variable and constants
//------------------------------------------------------------------------------


#ifndef __EXPREVAL_VALLIST_H
#define __EXPREVAL_VALLIST_H

// Includes
#include <string>
#include <vector>

// Part of expreval namespace
namespace ExprEval
    {
    // Value list item
    //--------------------------------------------------------------------------
    class ValueListItem
        {
        public:
            ValueListItem(const ::std::string &name, double def = 0.0, bool constant = false);
            ValueListItem(const ::std::string &name, double *ptr, double def = 0.0, bool constant = false);
            
            const ::std::string &GetName() const;
            bool IsConstant() const;
            
            double *GetAddress();
            void Reset();
            
        private:
            ::std::string m_name; // Name of value
            bool m_constant; // Value is constant
            
            double m_value; // Internal value (if ptr == 0)
            double *m_ptr; // Pointer to extern value if not 0
                
            double m_def; // Default value when reset
        };
        
    
        
        
    // Value list
    //--------------------------------------------------------------------------
    class ValueList
        {
        public:
            typedef ::std::vector<ValueListItem*>::size_type size_type;
                    
            ValueList();
            ~ValueList();
            
            // Add variable or constant to the list
            void Add(const ::std::string &name, double def = 0.0, bool constant = false);
            
            // Add an external variable or constant to the list
            void AddAddress(const ::std::string &name, double *ptr, double def = 0.0, bool constant = false);
            
            // Get the address of a variable or constant, internal or external
            double *GetAddress(const ::std::string &name) const;
            
            // Is the value constant
            bool IsConstant(const ::std::string &name) const;
            
            // Enumerate values
            size_type Count() const;
            void Item(size_type pos, ::std::string *name = 0, double *value = 0) const;
            
            // Initialize some default values (math constants)
            void AddDefaultValues();
            
            // Reset items to default values (constants are not changed)
            void Reset();
            
            // Free items and clear memory
            void Clear();
        
        private:
            ::std::vector<ValueListItem*> m_values;
        };
    }

#endif // __EXPREVAL_VALLIST_H

