#include "../include/RHF.h"

#define STRING_MAX 50

int main(int argc, char const *argv[])
{

    //FILES
    FILE * basisfile = NULL;
    FILE * inputfile = NULL;

    //loops
    int i,j;

    int length; // the total number of the orbitals -> the width/height of the matrices
    int el_num; // the total number of the electrons
    int el_num_counter; // counter for the loops concerning electrons

    el_num = 0;

    char basis[STRING_MAX]; // file dir for basis
    char input[STRING_MAX]; // file dir for input file
    char reader[STRING_MAX]; // reader for interpretation of input file

    double coord_temp; //help converting the cartisian coordinates from Angstrom unit into atomic unit

    strcpy(basis,"basis/");
    strcpy(input,"NULL");

    orbital * orbitals /*containing all orbitals*/, * orbital_temp;
    atomic_orbital * basis_HEAD; // containing all information from the basis
    atomic_orbital * atoms; // containing all information of the all atoms
    atomic_orbital * atoms_temp, * atoms_bk; // help copying the atoms from basis
    atomic_orbital * basis_scanner;


    gsl_matrix * overlap; //overlap matrix


    for(i=1;i<argc;i++)
    {
        strcpy(input,argv[i]);

        inputfile = fopen(input,"r");

        // interpret input file
        while(fscanf(inputfile,"%s",reader)!=EOF)
        {
            // read basis
            if(strcmp(reader,"&BASIS")==0)
            {
                // read the name of the basis
                fscanf(inputfile,"%s",reader);
                // concatenate it into basis dir
                strcat(basis,reader);
                // add postfix (to read)
                strcat(basis,".txt");
                // open file
                basisfile = fopen(basis,"r");
                //throw error if there is no basis
                if(basisfile==NULL)
                {
                    printf("HF_ERROR: BASIS NOT SUPPORTED!");
                    return 2;
                }
                //allocate memory for the basis_HEAD
                basis_HEAD = atomic_orbital_calloc();
                //read basis!
                basis_fscanf(basisfile,basis_HEAD);
                fclose(basisfile);
            }

            // read atoms' information
            if(strcmp(reader,"&COORD")==0)
            {
                // allocate -> set the head of the linked list
                atoms = atomic_orbital_calloc();

                // saving the location of the head
                atoms_bk = atoms;

                // read atoms
                while(fscanf(inputfile,"%s",reader)!=EOF)
                {
                    if(strcmp(reader,"&END_COORD")==0)
                    {
                        atomic_orbital_free(atoms);
                        atoms = atoms_bk;
                        atoms_temp->NEXT = NULL;
                        break;
                    }

                    else
                    {
                        // scanner for basis
                        basis_scanner = basis_HEAD;

                        // scan the basis
                        while(basis_scanner->NEXT!=NULL)
                        {
                            // whether the name matches
                            if(strcmp(basis_scanner->name,reader)==0)
                            {
                                atomic_orbital_single_cpy(atoms,basis_scanner);

                                break;
                            }

                            basis_scanner = basis_scanner->NEXT;
                        }

                        if(strcmp(basis_scanner->name,reader)==0)
                        {
                            // copy information from the basis
                            atomic_orbital_single_cpy(atoms,basis_scanner);
                            // add electrons
                            el_num += basis_scanner->N;
                            // copy coordinates & convert unit
                            for(i=0;i<3;i++)
                            {
                                fscanf(inputfile,"%lf",&coord_temp);
                                atoms->cartesian[i] = coord_temp / 0.529177210903;
                            }

                            // copy the coordinates of the atom to all of its orbitals
                            atomic_orbital_sync_coord(atoms);
                            // print the name of the atom to all its orbitals
                            atomic_orbital_name_print(atoms);

                            atoms->NEXT = atomic_orbital_calloc();

                            atoms_temp = atoms;

                            atoms = atoms->NEXT;                        
                        }

                        // no such an atom in the basis -> throw error
                        else
                        {
                            printf("OVERLAP_ERROR: atom %s is not defined in the basis.\n",reader);

                            return 3;
                        }
                        

                    }
                }
            }
        }

        // setting the head for the orbitals 
        orbitals = orbital_calloc((atoms->orbital_HEAD)->total);

        // save the location of the head
        orbital_temp = orbitals;
        // set the scanner of the atoms' list to its head
        atoms_temp = atoms;

        // copy all the orbitals from all atoms into `orbitals`
        while(atoms_temp->NEXT != NULL)
        {
            orbital_cpy(orbital_temp,atoms_temp->orbital_HEAD);
            while(orbital_temp->NEXT != NULL)
                orbital_temp = orbital_temp->NEXT;

            atoms_temp = atoms_temp->NEXT;
            orbital_temp->NEXT = orbital_calloc((atoms_temp->orbital_HEAD)->total);

            orbital_temp = orbital_temp->NEXT;
        }
        // DO the copying for the last one atom
        orbital_cpy(orbital_temp,atoms_temp->orbital_HEAD);

        // count the orbitals
        length = orbital_count(orbitals);

        overlap = gsl_matrix_calloc(length,length);

        orbital_S_matrix(overlap,orbitals);

        printf("Overlap matrix:\n");

        gsl_matrix_printf(overlap,length,length,"%14.6f");

        el_num_counter = 0;

        printf("\n\n");

        //print the labels of each atomic orbital in the order how they participate in the matrices
        printf("MO_LABEL:\n");
        printf("[");
        orbital_temp = orbitals;
        while(orbital_temp->NEXT != NULL)
        {
            printf(" %s ,",orbital_temp->label);
            orbital_temp = orbital_temp->NEXT;
        }
        printf(" %s ]\n\n",orbital_temp->label);


    }
    return 0;
}
