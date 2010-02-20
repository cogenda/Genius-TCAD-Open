// File:    parser.h
// Author:  Brian Vanderburg II
// Purpose: Parser object to help parse expression
//------------------------------------------------------------------------------


#ifndef __EXPREVAL_PARSER_H
#define __EXPREVAL_PARSER_H

// Includes
#include <vector>
#include <string>

// Part of expreval namespace
namespace ExprEval
    {
    // Forward declarations
    class Expression;
    class Node;
    
    // Token class
    //--------------------------------------------------------------------------
    class Token
        {
        public:
            
            // Token type
            enum TokenType
                {
                TypeUnknown = 0,
                
                // Basic types
                TypeOpenParenthesis, // '('
                TypeCloseParenthesis, // ')'
                TypeEqual, // '='
                TypePlus, // '+'
                TypeHyphen, // '-'
                TypeAsterisk, // '*'
                TypeForwardSlash, // '/'
                TypeHat, // '^'                
                TypeAmpersand, // '&'
                TypeAt, // @
                TypeComma, // ','                
                TypeSemicolon, // ';'
                TypeIdentifier, // name
                TypeValue // 0.2, .6, 2.1
                };            
        
            // Constructors
            Token(TokenType type, ::std::string::size_type start, ::std::string::size_type end);
            Token(const ::std::string &ident, ::std::string::size_type start, ::std::string::size_type end); // May throw if can not copy string
            Token(double value, ::std::string::size_type start, ::std::string::size_type end);
            
            // Information
            TokenType GetType() const;
            const ::std::string &GetIdentifier() const;
            double GetValue() const;
            
            ::std::string::size_type GetStart() const;
            ::std::string::size_type GetEnd() const;
            
        private:
            TokenType m_type; // Type of token
            ::std::string m_ident; // If type is Identifier
            double m_value; // If type is Value
            
            ::std::string::size_type m_start, m_end;
        };
        
    
        
    // Parser class
    //--------------------------------------------------------------------------
    class Parser
        {
        public:
            typedef ::std::vector<Token*>::size_type size_type;
        
            Parser(Expression *expr);
            ~Parser();
                        
            Node *Parse(const ::std::string &exstr);
            Node *ParseRegion(size_type start, size_type end);
            
            const Token &operator[] (size_type pos) const;
            
        private:

            void BuildTokens(const std::string &exstr);
        
            Expression *m_expr;
            ::std::vector<Token> m_tokens; // Store token and not pointer to token
        };
    }
    
#endif // __EXPREVAL_PARSER_H  

