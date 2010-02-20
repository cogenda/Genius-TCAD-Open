// File:    parser.cpp
// Author:  Brian Vanderburg II
// Purpose: Parser object to help parse expression
//------------------------------------------------------------------------------


// Includes
#include <new>
#include <memory>
#include <cstdlib>

#include "expr_parser.h"
#include "expr_node.h"
#include "expr_except.h"
#include "expr_funclist.h"
#include "expr.h"


using namespace std;
using namespace ExprEval;


// Token
//------------------------------------------------------------------------------

// Constructor
Token::Token(TokenType type, string::size_type start, string::size_type end) :
        m_type(type),
        m_start(start),
        m_end(end)
    {
    }
    
// Construct identifier token
Token::Token(const string &ident, string::size_type start, string::size_type end) : 
        m_type(Token::TypeIdentifier),
        m_ident(ident),
        m_start(start),
        m_end(end)
    {
    }
    
// Construct value token
Token::Token(double value, string::size_type start, string::size_type end) :
        m_type(Token::TypeValue),
        m_value(value),
        m_start(start),
        m_end(end)
    {
    }
    
// Get type
Token::TokenType Token::GetType() const
    {
    return m_type;
    }
    
// Get identifier
const string &Token::GetIdentifier() const
    {
    return m_ident;
    }
    
// Get value
double Token::GetValue() const
    {
    return m_value;
    }

// Get start
string::size_type Token::GetStart() const
    {
    return m_start;
    }
    
string::size_type Token::GetEnd() const
    {
    return m_end;
    }    

// Parser
//------------------------------------------------------------------------------

// Constructor
Parser::Parser(Expression *expr) : m_expr(expr)
    {
    if(expr == 0)
        throw(NullPointerException("Parser::Parser"));
    }
    
// Destructor
Parser::~Parser()
    {
    }
    
// Parse an expression string
Node *Parser::Parse(const string &exstr)
    {
    BuildTokens(exstr);
    
    // Make sure it is not still empty
    if(m_tokens.size() == 0)
        throw(EmptyExpressionException());
        
    // Parse the range
    return ParseRegion(0, m_tokens.size() - 1);
    }

// Parse a region of tokens
Node *Parser::ParseRegion(Parser::size_type start, Parser::size_type end)
    {
    size_type pos;
    size_type fgopen = (size_type)-1;
    size_type fgclose = (size_type)-1;
    size_type assignindex = (size_type)-1;
    size_type addsubindex = (size_type)-1;
    size_type muldivindex = (size_type)-1;
    size_type posnegindex = (size_type)-1;
    size_type expindex = (size_type)-1;
    bool multiexpr = false;
    int plevel = 0;
    
    // Check simple syntax
    if(start > end)
        throw(SyntaxException());
        
    // Scan through tokens
    for(pos = start; pos <= end; pos++)
        {
        switch(m_tokens[pos].GetType())
            {
            case Token::TypeOpenParenthesis:
                {
                plevel++;
                
                // Opening of first group?
                if(plevel == 1 && fgopen == (size_type)-1)
                    fgopen = pos;
                    
                break;
                };
                
            case Token::TypeCloseParenthesis:
                {
                plevel--;
                
                // First group closed?
                if(plevel == 0 && fgclose == (size_type)-1)
                    fgclose = pos;
                    
                if(plevel < 0)
                    {
                    UnmatchedParenthesisException e;
                    
                    e.SetStart(m_tokens[pos].GetStart());
                    e.SetEnd(m_tokens[pos].GetEnd());
                    throw(e);
                    }
                    
                break;
                }
                
            case Token::TypeEqual:
                {
                if(plevel == 0)
                    {
                    if(assignindex == (size_type)-1)
                        assignindex = pos;
                    }
                    
                break;
                }
                
            case Token::TypeAsterisk:
            case Token::TypeForwardSlash:
                {
                if(plevel == 0)
                    {
                    muldivindex = pos;
                    }
                    
                break;
                }
                
            case Token::TypeHat:
                {
                if(plevel == 0)
                    {
                    expindex = pos;
                    }
                
                break;
                }
                
            case Token::TypePlus:
            case Token::TypeHyphen:
                {
                if(plevel == 0)
                    {
                    if(pos == start)
                        {
                        // Positive or negative sign
                        if(posnegindex == (size_type)-1)
                            posnegindex = pos;
                        }
                    else
                        {
                        // What is before us
                        switch(m_tokens[pos - 1].GetType())
                            {
                            case Token::TypeEqual:
                            case Token::TypePlus:
                            case Token::TypeHyphen:
                            case Token::TypeAsterisk:
                            case Token::TypeForwardSlash:
                            case Token::TypeHat:
                                // After any of these, we are a positive/negative
                                if(posnegindex == (size_type)-1)
                                    posnegindex = pos;
                                break;
                                
                            default:
                                // After any other, we are addition/subtration
                                addsubindex = pos;
                                break;
                            }
                        }
                    }
                    
                break;
                }
                
            case Token::TypeSemicolon:
                {
                if(plevel == 0)
                    {
                    multiexpr = true;
                    }
                
                break;
                }
            }
        }
        
    // plevel should be 0
    if(plevel != 0)
        {
        UnmatchedParenthesisException e;
        
        e.SetStart(end);
        e.SetEnd(end);
        throw(e);
        }
        
    // Parse in certain order to maintain order of operators
    
    // Multi-expression first
    if(multiexpr)
        {
        auto_ptr<Node> n(new MultiNode(m_expr));
        n->Parse(*this, start, end);
        return n.release();
        }
    else if(assignindex != (size_type)-1)
        {
        // Assignment next
        auto_ptr<Node> n(new AssignNode(m_expr));
        n->Parse(*this, start, end, assignindex);
        return n.release();
        }
    else if(addsubindex != (size_type)-1)
        {
        // Addition/subtraction next
        if(m_tokens[addsubindex].GetType() == Token::TypePlus)
            {
            // Addition
            auto_ptr<Node> n(new AddNode(m_expr));
            n->Parse(*this, start, end, addsubindex);
            return n.release();
            }
        else
            {
            // Subtraction
            auto_ptr<Node> n(new SubtractNode(m_expr));
            n->Parse(*this, start, end, addsubindex);
            return n.release();
            }
        }
    else if(muldivindex != (size_type)-1)
        {
        // Multiplication/division next
        
        if(m_tokens[muldivindex].GetType() == Token::TypeAsterisk)
            {
            // Multiplication
            auto_ptr<Node> n(new MultiplyNode(m_expr));
            n->Parse(*this, start, end, muldivindex);
            return n.release();
            }
        else
            {
            // Division
            auto_ptr<Node> n(new DivideNode(m_expr));
            n->Parse(*this, start, end, muldivindex);
            return n.release();
            }
        }
    else if(posnegindex == start)
        {
        // Positive/negative next, must be at start and check before exponent
        if(m_tokens[posnegindex].GetType() == Token::TypePlus)
            {
            // Positive
            return ParseRegion(posnegindex + 1, end);
            }
        else
            {
            auto_ptr<Node> n(new NegateNode(m_expr));
            n->Parse(*this, start, end, posnegindex);
            return n.release();
            }
        }
    else if(expindex != (size_type)-1)
        {
        // Exponent
        auto_ptr<Node> n(new ExponentNode(m_expr));
        n->Parse(*this, start, end, expindex);
        return n.release();
        }
    else if(posnegindex != (size_type)-1)
        {
        // Check pos/neg again.  After testing for exponent, a pos/neg
        // at plevel 0 is syntax error
        SyntaxException e;
        
        e.SetStart(m_tokens[posnegindex].GetStart());
        e.SetEnd(m_tokens[posnegindex].GetEnd());
        throw(e);
        }
    else if(fgopen == start)
        {
        // Group parenthesis, make sure something in between them
        if(fgclose == end && fgclose > fgopen + 1)
            {
            return ParseRegion(fgopen + 1, fgclose - 1);
            }
        else
            {
            SyntaxException e;
            
            e.SetStart(m_tokens[fgopen].GetStart());
            if(fgclose == (size_type)-1)
                e.SetEnd(m_tokens[fgopen].GetEnd());
            else
                e.SetEnd(m_tokens[fgclose].GetEnd());
                
            throw(e);
            }
        }
    else if(fgopen == start + 1)
        {
        // Function
        if(fgclose == end)
            {
            // Find function list
            FunctionList *flist = m_expr->GetFunctionList();
            
            if(flist == 0)
                {
                NoFunctionListException e;
        
                e.SetStart(m_tokens[start].GetStart());
                e.SetEnd(m_tokens[start].GetEnd());
                throw(e);
                }
                
            // Get name
            string ident = m_tokens[start].GetIdentifier();
            
            // Create function node
            auto_ptr<FunctionNode> n(flist->Create(ident, m_expr));
            
            if(n.get())
                {
                n->Parse(*this, fgopen, fgclose);
                }
            else
                {
                NotFoundException e(ident);
        
                e.SetStart(m_tokens[start].GetStart());
                e.SetEnd(m_tokens[start].GetEnd());
                throw(e);                
                }
                
            return n.release();
            }
        else
            {
            SyntaxException e;
        
            e.SetStart(m_tokens[fgopen].GetStart());
            if(fgclose == (size_type)-1)
                e.SetEnd(m_tokens[fgopen].GetEnd());
            else
                e.SetEnd(m_tokens[fgclose].GetEnd());    
                
            throw(e);
            }
        }
    else if(start == end)
        {
        // Value, variable, or constant
        
        if(m_tokens[start].GetType() == Token::TypeIdentifier)
            {
            // Variable/constant
            auto_ptr<Node> n(new VariableNode(m_expr));
            n->Parse(*this, start, end);
            return n.release();    
            }
        else
            {
            // Value
            auto_ptr<Node> n(new ValueNode(m_expr));
            n->Parse(*this, start, end);
            return n.release();
            }
        }
    else
        {
        // Unknown, syntax
        SyntaxException e;
        
        e.SetStart(m_tokens[pos].GetStart());
        e.SetEnd(m_tokens[pos].GetEnd());
        
        throw(e);
        }
    }
    
// Get a token
const Token &Parser::operator[] (Parser::size_type pos) const
    {
    return m_tokens[pos];
    }
    
// Build tokens
void Parser::BuildTokens(const string &exstr)
    {
    m_tokens.clear();
    
    // Test zero-length expression
    if(exstr.length() == 0)
        {
        throw(EmptyExpressionException());
        }
        
    // Search through list
    string::size_type pos;
    bool comment = false;
    
    for(pos = 0; pos < exstr.length(); pos++)
        {
        // Take action based on character
        switch(exstr[pos])
            {
            // Comment
            case '#':
                {
                comment = true;
                    
                break;
                }
                
            // Newline ends comment
            case '\r':
            case '\n':
                {
                comment = false;
                
                break;
                }
                
            // Open parenthesis
            case '(':
                {
                if(comment == false)
                    {
                    m_tokens.push_back(Token(Token::TypeOpenParenthesis, pos, pos));
                    }
                break;
                }
                
            // Close parenthesis
            case ')':
                {
                if(comment == false)
                    {
                    m_tokens.push_back(Token(Token::TypeCloseParenthesis, pos, pos));
                    }
                break;
                }
                
            // Equal
            case '=':
                {
                if(comment == false)
                    {
                    m_tokens.push_back(Token(Token::TypeEqual, pos, pos));
                    }
                break;
                }
                
            // Plus
            case '+':
                {
                if(comment == false)
                    {
                    m_tokens.push_back(Token(Token::TypePlus, pos, pos));
                    }
                break;
                }
                
            // Hyphen
            case '-':
                {
                if(comment == false)
                    {
                    m_tokens.push_back(Token(Token::TypeHyphen, pos, pos));
                    }
                break;
                }
                
            // Asterisk
            case '*':
                {
                if(comment == false)
                    {
                    m_tokens.push_back(Token(Token::TypeAsterisk, pos, pos));
                    }
                break;
                }
                
            // Forward slash
            case '/':
                {
                if(comment == false)
                    {
                    m_tokens.push_back(Token(Token::TypeForwardSlash, pos, pos));
                    }
                break;
                }
                
            // Hat (exponent)
            case '^':
                {
                if(comment == false)
                    {
                    m_tokens.push_back(Token(Token::TypeHat, pos, pos));
                    }
                break;
                }
                
            // Ampersand
            case '&':
                {
                if(comment == false)
                    {
                    m_tokens.push_back(Token(Token::TypeAmpersand, pos, pos));
                    }
                break;
                }
                
            // At sign 
            case '@':
                {
                if(comment == false)
                    {
                    m_tokens.push_back(Token(Token::TypeAt, pos, pos));
                    }
                break;    
                }                
                
            // Comma
            case ',':
                {
                if(comment == false)
                    {
                    m_tokens.push_back(Token(Token::TypeComma, pos, pos));
                    }
                break;
                }
                
            // Semicolon
            case ';':
                {
                if(comment == false)
                    {
                    m_tokens.push_back(Token(Token::TypeSemicolon, pos, pos));
                    }
                break;
                }
                
            // None of the above, but it may be an identifier or value
            default:
                {
                if(comment == false)
                    {
                    // First, test for value
                    if(exstr[pos] == '.' || isdigit(exstr[pos]))
                        {
                        // We are a value
                        string::size_type start = pos;
                        
                        // Digits before period
                        while(isdigit(exstr[pos]))
                            pos++;
                            
                        // Period
                        if(exstr[pos] == '.')
                            pos++;
                            
                        // Digits after period
                        while(isdigit(exstr[pos]))
                            pos++;
                            
                        // Create token
                        string ident = exstr.substr(start, pos - start);
                        m_tokens.push_back(Token(atof(ident.c_str()), start, pos - 1));
                        
                        // Move pos back so pos++ will set it right
                        pos--;
                        }
                    else if(exstr[pos] == '_' || isalpha(exstr[pos]))
                        {
                        // We are an identifier
                        string::size_type start = pos;
                        bool foundname = true; // Found name part
                        
                        // Search for name, then period, etc
                        // An identifier can be multiple parts.  Each part
                        // is formed as an identifier and seperated by a period,
                        // An identifier can not end in a period
                        //
                        // color1.red : 1 identifier token
                        // color1. : color1 is identifier, . begins new token
                        // color1.1red : Not value (part 2 is not right)
                        while(foundname)
                            {
                            // Part before period
                            while(exstr[pos] == '_' || isalnum(exstr[pos]))
                                pos++;
                                
                            // Is there a period
                            if(exstr[pos] == '.')
                                {
                                pos++;
                                
                                // There is a period, look for the name again
                                if(exstr[pos] == '_' || isalpha(exstr[pos]))
                                    {
                                    foundname = true;
                                    }
                                else
                                    {
                                    // No name after period
                                    foundname = false;
                                    
                                    // Remove period from identifier
                                    pos--;
                                    }
                                }
                            else
                                {
                                // No period after name, so no new name
                                foundname = false;
                                }
                            }
                            
                        // Create token
                        m_tokens.push_back(Token(exstr.substr(start, pos - start), start, pos - 1));
                        
                        // Move pos back so pos++ will set it right
                        pos--;
                        }
                    else if(isspace(exstr[pos]))
                        {
                        // Do nothing, just ignore white space, but it still
                        // seperates tokens
                        }
                    else
                        {
                        // Unknown token
                        UnknownTokenException e;
                        e.SetStart(pos);
                        e.SetEnd(pos);
                        
                        throw(e);
                        }
                    }
                break;
                }
            }
        }
    }    
    
    
