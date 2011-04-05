// File:    node.cpp
// Author:  Brian Vanderburg II
// Purpose: Expression node
//------------------------------------------------------------------------------


// Includes
#include <new>
#include <memory>
#include <cmath>
#include <cerrno>

#include "expr_node.h"
#include "expr.h"
#include "expr_vallist.h"
#include "expr_funclist.h"
#include "expr_datalist.h"
#include "expr_except.h"

using namespace std;
using namespace ExprEval;

// Node
//------------------------------------------------------------------------------

// Constructor
Node::Node(Expression *expr) : m_expr(expr)
    {
    if(expr == 0)
        throw(NullPointerException("Node::Node"));
    }

// Destructor
Node::~Node()
    {
    }

// Evaluate
double Node::Evaluate()
    {
    m_expr->TestAbort();

    return DoEvaluate();
    }

// Function node
//------------------------------------------------------------------------------

// Constructor
FunctionNode::FunctionNode(Expression *expr) : Node(expr),
        m_argMin(0), m_argMax(0), m_refMin(0), m_refMax(0), m_dataMin(0), m_dataMax(0)
    {
    }

// Destructor
FunctionNode::~FunctionNode()
    {
    // Delete child nodes
    vector<Node*>::size_type pos;

    for(pos = 0; pos < m_nodes.size(); pos++)
        {
        delete m_nodes[pos];
        }
    }

// Get name
string FunctionNode::GetName() const
    {
    return m_factory->GetName();
    }

// Set argument count
void FunctionNode::SetArgumentCount(long argMin, long argMax, long refMin, long refMax,
        long dataMin, long dataMax)
    {
    m_argMin = argMin;
    m_argMax = argMax;
    m_refMin = refMin;
    m_refMax = refMax;
    m_dataMin = dataMin;
    m_dataMax = dataMax;
    }

// Parse expression
void FunctionNode::Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
        Parser::size_type v1)
    {
    Parser::size_type pos, last;
    int plevel = 0;

    // Sanity check (start/end are function parenthesis)
    if(start >= end)
        throw(SyntaxException());

    // Look
    last = start + 1;
    for(pos = start + 1; pos <= end && start + 1 != end; pos++)
        {
        switch(parser[pos].GetType())
            {
            case Token::TypeOpenParenthesis:
                {
                plevel++;
                break;
                };

            case Token::TypeComma:
            case Token::TypeCloseParenthesis:
                {
                // Handle Close parenthesis for all but the last one at the end
                if(parser[pos].GetType() == Token::TypeCloseParenthesis && pos != end)
                    {
                    plevel--;

                    if(plevel < 0)
                        {
                        UnmatchedParenthesisException e;

                        e.SetStart(parser[pos].GetStart());
                        e.SetEnd(parser[pos].GetEnd());

                        throw(e);
                        }

                    break;
                    }

                // Handle comma, or if it was the ending parenthesis treat it like comma
                if(plevel == 0)
                    {
                    if(pos > last)
                        {
                        // reference parameter?
                        if(parser[last].GetType() == Token::TypeAmpersand)
                            {
                            // Reference parameter, check position and type of next parameter
                            if(last == pos - 2 && parser[last + 1].GetType() == Token::TypeIdentifier)
                                {
                                // Get value list
                                ValueList *vlist = m_expr->GetValueList();
                                if(vlist == 0)
                                    {
                                    NoValueListException e;

                                    e.SetStart(parser[last + 1].GetStart());
                                    e.SetEnd(parser[last + 1].GetEnd());
                                    throw(e);
                                    }

                                // Get name
                                string ident = parser[last + 1].GetIdentifier();

                                // Make sure it is not a constant
                                if(vlist->IsConstant(ident))
                                    {
                                    ConstantReferenceException e(ident);

                                    e.SetStart(parser[last + 1].GetStart());
                                    e.SetEnd(parser[last + 1].GetEnd());
                                    throw(e);
                                    }

                                // Get address
                                double *vaddr = vlist->GetAddress(ident);
                                if(vaddr == 0)
                                    {
                                    // Try to add it and get again
                                    vlist->Add(ident);
                                    vaddr = vlist->GetAddress(ident);
                                    }

                                if(vaddr == 0)
                                    {
                                    NotFoundException e(ident);

                                    e.SetStart(parser[last + 1].GetStart());
                                    e.SetEnd(parser[last + 1].GetEnd());

                                    throw(e);
                                    }

                                // Add it
                                m_refs.push_back(vaddr);
                                }
                            else
                                {
                                SyntaxException e;

                                e.SetStart(parser[last].GetStart());
                                e.SetEnd(parser[pos].GetEnd());
                                throw(e);
                                }
                            } // TypeAmpersand
                        else if(parser[last].GetType() == Token::TypeAt)
                            {
                            // Data parameter, check position and type of next parameter
                            if(last == pos - 2 && parser[last + 1].GetType() == Token::TypeIdentifier)
                                {
                                // Get data list
                                DataList *dlist = m_expr->GetDataList();
                                if(dlist == 0)
                                    {
                                    NoDataListException e;

                                    e.SetStart(parser[last + 1].GetStart());
                                    e.SetEnd(parser[last + 1].GetEnd());
                                    throw(e);
                                    }

                                // Get name
                                string ident = parser[last + 1].GetIdentifier();

                                // Get address
                                DataEntry *daddr = dlist->GetDataEntry(ident);
                                if(daddr == 0)
                                    {
                                    // Try to add it and get again
                                    dlist->Add(ident);
                                    daddr = dlist->GetDataEntry(ident);
                                    }

                                if(daddr == 0)
                                    {
                                    NotFoundException e(ident);

                                    e.SetStart(parser[last + 1].GetStart());
                                    e.SetEnd(parser[last + 1].GetEnd());

                                    throw(e);
                                    }

                                // Add it
                                m_data.push_back(daddr);
                                }
                            else
                                {
                                SyntaxException e;

                                e.SetStart(parser[last].GetStart());
                                e.SetEnd(parser[pos].GetEnd());
                                throw(e);
                                }
                            } // TypeAt
                        else
                            {
                            // Create node
                            auto_ptr<Node> n(parser.ParseRegion(last, pos - 1));
                            m_nodes.push_back(n.get());
                            n.release();
                            }
                        }
                    else
                        {
                        SyntaxException e;

                        e.SetStart(parser[pos].GetStart());
                        e.SetEnd(parser[pos].GetEnd());
                        throw(e);
                        }

                    last = pos + 1;
                    }

                break;
                }
            }
        }

    // plevel should be zero
    if(plevel != 0)
        {
        UnmatchedParenthesisException e;

        e.SetStart(parser[end].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }


    // Check argument count
    if(m_argMin != -1 && m_nodes.size() < (vector<Node*>::size_type)m_argMin)
        {
        InvalidArgumentCountException e(GetName());

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    if(m_argMax != -1 && m_nodes.size() > (vector<Node*>::size_type)m_argMax)
        {
        InvalidArgumentCountException e(GetName());

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    if(m_refMin != -1 && m_refs.size() < (vector<double*>::size_type)m_refMin)
        {
        InvalidArgumentCountException e(GetName());

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    if(m_refMax != -1 && m_refs.size() > (vector<double*>::size_type)m_refMax)
        {
        InvalidArgumentCountException e(GetName());

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    if(m_dataMin != -1 && m_data.size() < (vector<DataEntry*>::size_type)m_dataMin)
        {
        InvalidArgumentCountException e(GetName());

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    if(m_dataMax != -1 && m_data.size() > (vector<DataEntry*>::size_type)m_dataMax)
        {
        InvalidArgumentCountException e(GetName());

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }
    }

// Multi node
//------------------------------------------------------------------------------

// Constructor
MultiNode::MultiNode(Expression *expr) : Node(expr)
    {
    }

// Destructor
MultiNode::~MultiNode()
    {
    // Free child nodes
    vector<Node*>::size_type pos;

    for(pos = 0; pos < m_nodes.size(); pos++)
        {
        delete m_nodes[pos];
        }
    }

// Evaluate
double MultiNode::DoEvaluate()
    {
    vector<Node*>::size_type pos;
    double result = 0.0;

    for(pos = 0; pos < m_nodes.size(); pos++)
        {
        result = m_nodes[pos]->Evaluate();
        }

    return result;
    }

// Parse
void MultiNode::Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
        Parser::size_type v1)
    {
    Parser::size_type pos, last;
    int plevel = 0;

    // Sanity check
    if(start >= end)
        throw(SyntaxException());

    // Look
    last = start;
    for(pos = start; pos <= end; pos++)
        {
        switch(parser[pos].GetType())
            {
            case Token::TypeOpenParenthesis:
                {
                plevel++;
                break;
                };

            case Token::TypeCloseParenthesis:
                {
                plevel--;

                if(plevel < 0)
                    {
                    UnmatchedParenthesisException e;

                    e.SetStart(parser[pos].GetStart());
                    e.SetEnd(parser[pos].GetEnd());
                    throw(e);
                    }

                break;
                }

            case Token::TypeSemicolon:
                {
                if(plevel == 0)
                    {
                    if(pos > last)
                        {
                        // Everything from last to pos - 1
                        auto_ptr<Node> n(parser.ParseRegion(last, pos - 1));
                        m_nodes.push_back(n.get());
                        n.release();
                        }
                    else
                        {
                        SyntaxException e;

                        e.SetStart(parser[last].GetStart());
                        e.SetEnd(parser[pos].GetEnd());
                        throw(e);
                        }

                    last = pos + 1;
                    }

                break;
                }
            }
        }

    // plevel should be zero
    if(plevel != 0)
        {
        UnmatchedParenthesisException e;

        e.SetStart(parser[pos].GetStart());
        e.SetEnd(parser[pos].GetEnd());
        throw(e);
        }

    // If the end was not a semicolon, test it as well
    if(last < end + 1)
        {
        auto_ptr<Node> n(parser.ParseRegion(last, end));
        m_nodes.push_back(n.get());
        n.release();
        }
    }

// Assign node
//------------------------------------------------------------------------------

// Constructor
AssignNode::AssignNode(Expression *expr) : Node(expr), m_var(0), m_rhs(0)
    {
    }

// Destructor
AssignNode::~AssignNode()
    {
    // Free child node
    delete m_rhs;
    }

// Evaluate
double AssignNode::DoEvaluate()
    {
    return (*m_var = m_rhs->Evaluate());
    }

// Parse
void AssignNode::Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
        Parser::size_type v1)
    {
    // Check some basic syntax
    if(v1 != start + 1 || v1 >= end)
        {
        SyntaxException e;

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    // Get the value list
    ValueList *vlist = m_expr->GetValueList();
    if(vlist == 0)
        {
        NoValueListException e;

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[start].GetEnd());
        throw(e);
        }

    // Determine variable name
    string ident = parser[start].GetIdentifier();

    // Make sure it is not a constant
    if(vlist->IsConstant(ident))
        {
        ConstantAssignException e(ident);

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[v1].GetEnd());
        throw(e);
        }

    // Get address
    double *vaddr = vlist->GetAddress(ident);

    if(vaddr == 0)
        {
        // If it does not already exist, try to create it
        vlist->Add(ident);
        vaddr = vlist->GetAddress(ident);
        }

    if(vaddr == 0)
        {
        NotFoundException e(ident);

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[start].GetEnd());
        throw(e);
        }

    // Parse the node (will throw if it can not parse)
    auto_ptr<Node> n(parser.ParseRegion(v1 + 1, end));

    // Set data
    m_var = vaddr;
    m_rhs = n.release();
    }


// Add node
//------------------------------------------------------------------------------

// Constructor
AddNode::AddNode(Expression *expr) : Node(expr), m_lhs(0), m_rhs(0)
    {
    }

// Destructor
AddNode::~AddNode()
    {
    // Free child nodes
    delete m_lhs;
    delete m_rhs;
    }

// Evaluate
double AddNode::DoEvaluate()
    {
    return m_lhs->Evaluate() + m_rhs->Evaluate();
    }

// Parse
void AddNode::Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
        Parser::size_type v1)
    {
    // Check some basic syntax
    if(v1 <= start || v1 >= end)
        {
        SyntaxException e;

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    // Parse sides
    auto_ptr<Node> left(parser.ParseRegion(start, v1 - 1));
    auto_ptr<Node> right(parser.ParseRegion(v1 + 1, end));

    m_lhs = left.release();
    m_rhs = right.release();
    }

// Subtract node
//------------------------------------------------------------------------------

// Constructor
SubtractNode::SubtractNode(Expression *expr) : Node(expr), m_lhs(0), m_rhs(0)
    {
    }

// Destructor
SubtractNode::~SubtractNode()
    {
    // Free child nodes
    delete m_lhs;
    delete m_rhs;
    }

// Evaluate
double SubtractNode::DoEvaluate()
    {
    return m_lhs->Evaluate() - m_rhs->Evaluate();
    }

// Parse
void SubtractNode::Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
        Parser::size_type v1)
    {
    // Check some basic syntax
    if(v1 <= start || v1 >= end)
        {
        SyntaxException e;

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    // Parse sides
    auto_ptr<Node> left(parser.ParseRegion(start, v1 - 1));
    auto_ptr<Node> right(parser.ParseRegion(v1 + 1, end));

    m_lhs = left.release();
    m_rhs = right.release();
    }

// Multiply node
//------------------------------------------------------------------------------

// Constructor
MultiplyNode::MultiplyNode(Expression *expr) : Node(expr), m_lhs(0), m_rhs(0)
    {
    }

// Destructor
MultiplyNode::~MultiplyNode()
    {
    // Free child nodes
    delete m_lhs;
    delete m_rhs;
    }

// Evaluate
double MultiplyNode::DoEvaluate()
    {
    return m_lhs->Evaluate() * m_rhs->Evaluate();
    }

// Parse
void MultiplyNode::Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
        Parser::size_type v1)
    {
    // Check some basic syntax
    if(v1 <= start || v1 >= end)
        {
        SyntaxException e;

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    // Parse sides
    auto_ptr<Node> left(parser.ParseRegion(start, v1 - 1));
    auto_ptr<Node> right(parser.ParseRegion(v1 + 1, end));

    m_lhs = left.release();
    m_rhs = right.release();
    }

// Divide node
//------------------------------------------------------------------------------

// Constructor
DivideNode::DivideNode(Expression *expr) : Node(expr), m_lhs(0), m_rhs(0)
    {
    }

// Destructor
DivideNode::~DivideNode()
    {
    // Free child nodes
    delete m_lhs;
    delete m_rhs;
    }

// Evaluate
double DivideNode::DoEvaluate()
    {
    double r2 = m_rhs->Evaluate();

    if(r2 != 0.0)
        {
        return m_lhs->Evaluate() / r2;
        }
    else
        {
        throw(DivideByZeroException());
        }
    }

// Parse
void DivideNode::Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
        Parser::size_type v1)
    {
    // Check some basic syntax
    if(v1 <= start || v1 >= end)
        {
        SyntaxException e;

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    // Parse sides
    auto_ptr<Node> left(parser.ParseRegion(start, v1 - 1));
    auto_ptr<Node> right(parser.ParseRegion(v1 + 1, end));

    m_lhs = left.release();
    m_rhs = right.release();
    }

// Negate node
//------------------------------------------------------------------------------

// Constructor
NegateNode::NegateNode(Expression *expr) : Node(expr), m_rhs(0)
    {
    }

// Destructor
NegateNode::~NegateNode()
    {
    // Free child nodes
    delete m_rhs;
    }

// Evaluate
double NegateNode::DoEvaluate()
    {
    return -(m_rhs->Evaluate());
    }

// Parse
void NegateNode::Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
        Parser::size_type v1)
    {
    // Check some basic syntax
    if(start != v1 || v1 >= end)
        {
        SyntaxException e;

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    // Parse sides
    auto_ptr<Node> right(parser.ParseRegion(v1 + 1, end));

    m_rhs = right.release();
    }

// Exponent node
//------------------------------------------------------------------------------

// Constructor
ExponentNode::ExponentNode(Expression *expr) : Node(expr), m_lhs(0), m_rhs(0)
    {
    }

// Destructor
ExponentNode::~ExponentNode()
    {
    // Free child nodes
    delete m_lhs;
    delete m_rhs;
    }

// Evaluate
double ExponentNode::DoEvaluate()
    {
    errno = 0;

    double result =  pow(m_lhs->Evaluate(), m_rhs->Evaluate());

    if(errno)
        throw(MathException("^"));

    return result;
    }

// Parse
void ExponentNode::Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
        Parser::size_type v1)
    {
    // Check some basic syntax
    if(v1 <= start || v1 >= end)
        {
        SyntaxException e;

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    // Parse sides
    auto_ptr<Node> left(parser.ParseRegion(start, v1 - 1));
    auto_ptr<Node> right(parser.ParseRegion(v1 + 1, end));

    m_lhs = left.release();
    m_rhs = right.release();
    }

// Variable node
//------------------------------------------------------------------------------

// Constructor
VariableNode::VariableNode(Expression *expr) : Node(expr), m_var(0)
    {
    }

// Destructor
VariableNode::~VariableNode()
    {
    }

// Evaluate
double VariableNode::DoEvaluate()
    {
    return *m_var;
    }

// Parse
void VariableNode::Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
        Parser::size_type v1)
    {
    // Check some basic syntax
    if(start != end)
        {
        SyntaxException e;

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    // Get value list
    ValueList *vlist = m_expr->GetValueList();
    if(vlist == 0)
        {
        NoValueListException e;

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[start].GetEnd());
        throw(e);
        }

    // Get name
    string ident = parser[start].GetIdentifier();

    // Get address
    double *vaddr = vlist->GetAddress(ident);

    if(vaddr == 0)
        {
        // If it does not already exist, try to create it
        vlist->Add(ident);
        vaddr = vlist->GetAddress(ident);
        }

    if(vaddr == 0)
        {
        NotFoundException e(ident);

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[start].GetEnd());
        throw(e);
        }

    // Set information
    m_var = vaddr;
    }

// Value node
//------------------------------------------------------------------------------

// Constructor
ValueNode::ValueNode(Expression *expr) : Node(expr), m_val(0)
    {
    }

// Destructor
ValueNode::~ValueNode()
    {
    }

// Evaluate
double ValueNode::DoEvaluate()
    {
    return m_val;
    }

// Parse
void ValueNode::Parse(Parser &parser, Parser::size_type start, Parser::size_type end,
        Parser::size_type v1)
    {
    // Check basic syntax
    if(start != end)
        {
        SyntaxException e;

        e.SetStart(parser[start].GetStart());
        e.SetEnd(parser[end].GetEnd());
        throw(e);
        }

    // Set info
    m_val = parser[start].GetValue();
    }



