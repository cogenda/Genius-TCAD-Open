// File:    except.cc
// Author:  Brian Vanderburg II
// Purpose: ExprEval exceptions
//------------------------------------------------------------------------------


// Includes
#include "expr_except.h"

using namespace std;
using namespace ExprEval;


// Default/unknown ExprEval exception
//------------------------------------------------------------------------------
Exception::Exception() :
        m_start((string::size_type)-1),
        m_end((string::size_type)-1)
    {
    m_type = Type_Exception;
    }
    
Exception::~Exception() throw()
    {
    }    

Exception::Type Exception::GetType() const
    {
    return m_type;
    }
    
const string &Exception::GetValue() const
    {
    return m_value;
    }    
    
void Exception::SetStart(string::size_type start)
    {
    m_start = start;
    }
    
void Exception::SetEnd(string::size_type end)
    {
    m_end = end;
    }
    
string::size_type Exception::GetStart() const
    {
    return m_start;
    }
    
string::size_type Exception::GetEnd() const
    {
    return m_end;
    }                
    
// Not found exception
//------------------------------------------------------------------------------
NotFoundException::NotFoundException(const string &name)
    {
    m_type = Type_NotFoundException;
    m_value = name;
    }
    
// Already exists exception
//------------------------------------------------------------------------------
AlreadyExistsException::AlreadyExistsException(const string &name)
    {
    m_type = Type_AlreadyExistsException;
    m_value = name;
    }
    
// Null pointer exception
//------------------------------------------------------------------------------
NullPointerException::NullPointerException(const string &method)
    {
    m_type = Type_NullPointerException;
    m_value = method;
    }
    
// Math error exception
//------------------------------------------------------------------------------
MathException::MathException(const string &function)
    {
    m_type = Type_MathException;
    m_value = function;
    }
    
// Divide by zero
//------------------------------------------------------------------------------
DivideByZeroException::DivideByZeroException()
    {
    m_type = Type_DivideByZeroException;
    }
    
// No value list found during parsing
//------------------------------------------------------------------------------
NoValueListException::NoValueListException()
    {
    m_type = Type_NoValueListException;
    }
    
// No function list found during parsing
//------------------------------------------------------------------------------
NoFunctionListException::NoFunctionListException()
    {
    m_type = Type_NoFunctionListException;
    }

// No data list found during parsing
//------------------------------------------------------------------------------
NoDataListException::NoDataListException()
    {
    m_type = Type_NoDataListException;
    }    
    
// Expression abort
//------------------------------------------------------------------------------
AbortException::AbortException()
    {
    m_type = Type_AbortException;
    }
    
// Empty expression
//------------------------------------------------------------------------------
EmptyExpressionException::EmptyExpressionException()
    {
    m_type = Type_EmptyExpressionException;
    }
    
// Unknown token found during parsing
//------------------------------------------------------------------------------
UnknownTokenException::UnknownTokenException()
    {
    m_type = Type_UnknownTokenException;
    }
    
// Invalid argument count
//------------------------------------------------------------------------------
InvalidArgumentCountException::InvalidArgumentCountException(const string &function)
    {
    m_type = Type_InvalidArgumentCountException;
    m_value = function;
    }
    
// Assign to constant
//------------------------------------------------------------------------------
ConstantAssignException::ConstantAssignException(const string &value)
    {
    m_type = Type_ConstantAssignException;
    m_value = value;
    }
    
// Pass constant by reference
//------------------------------------------------------------------------------
ConstantReferenceException::ConstantReferenceException(const string &value)
    {
    m_type = Type_ConstantReferenceException;
    m_value = value;
    }
    
// Syntax error exception
//------------------------------------------------------------------------------
SyntaxException::SyntaxException()
    {
    m_type = Type_SyntaxException;
    }
    
// Unmatched parenthesis
//------------------------------------------------------------------------------
UnmatchedParenthesisException::UnmatchedParenthesisException()
    {
    m_type = Type_UnmatchedParenthesisException;
    }

// Invalid data
//------------------------------------------------------------------------------
InvalidDataException::InvalidDataException(const string &function)
    {
    m_type = Type_InvalidDataException;
    m_value = function;
    }
