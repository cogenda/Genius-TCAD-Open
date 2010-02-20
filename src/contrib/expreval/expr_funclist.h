// File:    funclist.h
// Author:  Brian Vanderburg II
// Purpose: Function list
//------------------------------------------------------------------------------


#ifndef __EXPREVAL_FUNCLIST_H
#define __EXPREVAL_FUNCLIST_H

// Includes
#include <string>
#include <vector>

// Part of expreval namespace
namespace ExprEval
    {
    // Forward declarations
    class FunctionNode;
    class Expression;
    
    // Function factory
    //--------------------------------------------------------------------------
    class FunctionFactory
        {
        public:
            FunctionFactory();
            virtual ~FunctionFactory();
            
            virtual ::std::string GetName() const = 0;
            virtual FunctionNode *DoCreate(Expression *expr) = 0;
            
            FunctionNode *Create(Expression *expr);
        };
        
    // Function list
    //--------------------------------------------------------------------------
    class FunctionList
        {
        public:
            typedef ::std::vector<FunctionFactory*>::size_type size_type;
                    
            FunctionList();
            ~FunctionList();
            
            // Add a function factory to the list
            void Add(FunctionFactory *factory);
            
            // Create a node for a function
            FunctionNode *Create(const ::std::string &name, Expression *expr);
            
            // Initialize default functions
            void AddDefaultFunctions();
            
            // Free items and clear memory
            void Clear();
        
        private:
            ::std::vector<FunctionFactory*> m_functions;
        };
        
    }

#endif // __EXPREVAL_FUNCLIST_H

