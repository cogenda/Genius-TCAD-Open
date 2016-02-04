#ifndef __ngspice_interface_h__
#define __ngspice_interface_h__

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * the circuit node structure
   */
  typedef struct sCKTnode
  {
    char *name;
    int  type;

#define SP_VOLTAGE 3
#define SP_CURRENT 4
#define NODE_VOLTAGE SP_VOLTAGE
#define NODE_CURRENT SP_CURRENT

    int number;			/* Number of the node */
    double ic;			/* Value of the initial condition */
    double nodeset;		/* Value of the .nodeset option */
    double *ptr;		/* ??? */
    struct sCKTnode *next;	/* pointer to the next node */
    unsigned int icGiven;	/* FLAG ic given */
    unsigned int nsGiven;	/* FLAG nodeset given */
  }
  CKTnode;


#ifdef SPICE_STATIC_INTERFACE
  /**
   * init the netlist from spice input file
   */
  extern int init_netlist(char * ckt_file);


  /**
   * init circuit internal data structure
   */
  extern int init_ckt();


  /**
   * @return the dof number in circuit
   * please note that node 0 is not considered in matrix assemble
   */
  extern unsigned int n_dofs();


  /**
   * @return the node number in circuit
   * including voltage node and current branch node
   */
  extern unsigned int n_nodes();


  /**
   * get nodal information
   */
  extern void get_nodal_info(CKTnode **nodes);


  /**
   * get the matrix
   */
  extern MatrixPtr get_matrix();

  /**
   * get the address of matrix entry (row,col)
   */
  extern double * get_matrix_entry(int row, int col);

#endif


#ifdef __cplusplus
}
#endif

#endif
