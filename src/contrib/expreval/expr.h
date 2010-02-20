// File:    expr.h
// Author:  Brian Vanderburg II
// Purpose: Expression object
//------------------------------------------------------------------------------


#ifndef __EXPREVAL_EXPR_H
#define __EXPREVAL_EXPR_H

// Includes
#include <string>

// Part of expreval namespace
namespace ExprEval
    {
    // Forward declarations
    class ValueList;
    class FunctionList;
    class DataList;
    class Node;
    
    // Expression class
    //--------------------------------------------------------------------------
    class Expression
        {
        public:
            Expression();
            virtual ~Expression();
            
            // Variable list
            void SetValueList(ValueList *vlist);
            ValueList *GetValueList() const;
            
            // Function list
            void SetFunctionList(FunctionList *flist);
            FunctionList *GetFunctionList() const;
            
            // Data list
            void SetDataList(DataList *dlist);
            DataList *GetDataList() const;
            
            // Abort control
            virtual bool DoTestAbort();
            void TestAbort(bool force = false);
            void SetTestAbortCount(unsigned long count);
            
            // Parse an expression
            void Parse(const ::std::string &exstr);
            
            // Clear an expression
            void Clear();
            
            // Evaluate expression
            double Evaluate();
            
        protected:
            ValueList *m_vlist;
            FunctionList *m_flist;
            DataList *m_dlist;
            Node *m_expr;
            unsigned long m_abortcount;
            unsigned long m_abortreset;
        };

       
        
    } // namespace ExprEval
    
#endif // __EXPREVAL_EXPR_H  

