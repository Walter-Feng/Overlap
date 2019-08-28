#include "../include/integral.h"

//Calculate the Overlap Integral of two gaussian function
double SIntegral(double ra[3], double rb[3], int ax, int ay, int az, int bx, int by, int bz, double alpha,double beta)
{

    // if one of the angular number goes below zero, it means it will not have contribution - the same as giving derivation to a constant
    if(ax<0||ay<0||az<0||bx<0||by<0||bz<0) return 0;

    //Provide recurrence relation
    else if(ax > 0) return ((alpha*ra[0]+beta*rb[0])/(alpha + beta)-ra[0])*SIntegral(ra,rb,ax-1,ay,az,bx,by,bz,alpha,beta) + (ax-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax-2,ay,az,bx,by,bz,alpha,beta) + bx/2.0/(alpha+beta) * SIntegral(ra,rb,ax-1,ay,az,bx-1,by,bz,alpha,beta);

    else if(ay > 0) return ((alpha*ra[1]+beta*rb[1])/(alpha + beta)-ra[1])*SIntegral(ra,rb,ax,ay-1,az,bx,by,bz,alpha,beta) + (ay-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay-2,az,bx,by,bz,alpha,beta) + by/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay-1,az,bx,by-1,bz,alpha,beta);

    else if(az > 0) return ((alpha*ra[2]+beta*rb[2])/(alpha + beta)-ra[2])*SIntegral(ra,rb,ax,ay,az-1,bx,by,bz,alpha,beta) + (az-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az-2,bx,by,bz,alpha,beta) + bz/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az-1,bx,by,bz-1,alpha,beta);

    else if(bx > 0) return ((alpha*ra[0]+beta*rb[0])/(alpha + beta)-rb[0])*SIntegral(ra,rb,ax,ay,az,bx-1,by,bz,alpha,beta) + ax/2.0/(alpha+beta) * SIntegral(ra,rb,ax-1,ay,az,bx-1,by,bz,alpha,beta) + (bx-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx-2,by,bz,alpha,beta);

    else if(by > 0) return ((alpha*ra[1]+beta*rb[1])/(alpha + beta)-rb[1])*SIntegral(ra,rb,ax,ay,az,bx,by-1,bz,alpha,beta) + ay/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay-1,az,bx,by-1,bz,alpha,beta) + (by-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx,by-2,bz,alpha,beta);

    else if(bz > 0) return ((alpha*ra[2]+beta*rb[2])/(alpha + beta)-rb[2])*SIntegral(ra,rb,ax,ay,az,bx,by,bz-1,alpha,beta) + az/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az-1,bx,by,bz-1,alpha,beta) + (bz-1)/2.0/(alpha+beta) * SIntegral(ra,rb,ax,ay,az,bx,by,bz-2,alpha,beta) ;

    //giving the starting point
    else return sqrt(M_PI/(alpha + beta)) * M_PI/(alpha + beta) * exp(- alpha * beta /(alpha + beta) * (pow(ra[0]-rb[0],2)+pow(ra[1]-rb[1],2)+pow(ra[2]-rb[2],2)));
}


//allocate memory for struct gaussian_chain
gaussian_chain * gaussian_chain_calloc()
{
    gaussian_chain * temp;

    int i;

    temp = new gaussian_chain;

    for(i=0;i<3;i++)
    {
        temp->R[i] = 0;
        temp->a[i] = 0;
    }
    temp->exponent = 0;
    temp->coefficient = 0;

    temp->NEXT = NULL;
    return temp;
}

//free memory for struct gaussian_chain (not suitable for a single gaussian function)
void gaussian_chain_free(gaussian_chain * HEAD)
{
    gaussian_chain * temp1, *temp2;

    temp1 = HEAD;

    while(temp1->NEXT != NULL)
    {
        temp2 = temp1;
        temp1 = temp1->NEXT;
        delete temp2;
    }
    delete temp1;
}

//transform the struct orbital to struct gaussian_chain
void single_electron_transform(gaussian_chain * HEAD, orbital * a)
{
    gaussian_chain * temp, * bk;

    temp = HEAD;

    int i,j,q;

    for(i=0;i<a->length;i++)
    {
        for(j=0;j<a->total;j++)
        {
            temp->coefficient = a->A[i].coef * *(a->coefficients + j) * normalize(*(a->exponents + j),a->A[i].a[0],a->A[i].a[1],a->A[i].a[2]);
            temp->exponent = *(a->exponents + j);
            for(q=0;q<3;q++)
            {
                temp->R[q] = a->cartesian[q];
                temp->a[q] = a->A[i].a[q];
            }
            temp->NEXT = gaussian_chain_calloc();
            bk = temp;
            temp = temp->NEXT;
        }
    }
    delete temp;
    bk->NEXT = NULL;
}

//enabling overlap integrals for gaussian_chain format
double gaussian_chain_SIntegral(gaussian_chain * a, gaussian_chain * b)
{
    return a->coefficient * b->coefficient *SIntegral(a->R,b->R,a->a[0],a->a[1],a->a[2],b->a[0],b->a[1],b->a[2],a->exponent,b->exponent);
}

//enabling overlap integrals for gaussian_chain format, calculating the whole chain
double gaussian_chain_full_SIntegral(gaussian_chain * a_HEAD, gaussian_chain * b_HEAD)
{
    double result = 0;
    gaussian_chain * a_temp, * b_temp;

    a_temp = a_HEAD;
    b_temp = b_HEAD;

    while(a_temp->NEXT != NULL)
    {
        b_temp = b_HEAD;
        while(b_temp->NEXT != NULL)
        {
            result += gaussian_chain_SIntegral(a_temp,b_temp);

            b_temp = b_temp->NEXT;
        }
        result += gaussian_chain_SIntegral(a_temp,b_temp);
        a_temp = a_temp->NEXT;
    }

    b_temp = b_HEAD;
    while(b_temp->NEXT != NULL)
    {
        result += gaussian_chain_SIntegral(a_temp,b_temp);

        b_temp = b_temp->NEXT;
    }

    result += gaussian_chain_SIntegral(a_temp,b_temp);

    return result;
}

//enabling ovelap integrals for orbital format
double orbital_SIntegral(orbital * a, orbital * b)
{
    double result;

    result = 0;
    gaussian_chain * a_head, * b_head, * a_temp, * b_temp;
    a_head = gaussian_chain_calloc();
    b_head = gaussian_chain_calloc();

    single_electron_transform(a_head,a);
    single_electron_transform(b_head,b);

    a_temp = a_head;

    while(a_temp->NEXT != NULL)
    {
        b_temp = b_head;
        while(b_temp->NEXT != NULL)
        {
            result += gaussian_chain_SIntegral(a_temp,b_temp);

            b_temp = b_temp->NEXT;
        }
        result += gaussian_chain_SIntegral(a_temp,b_temp);
        a_temp = a_temp->NEXT;
    }

    b_temp = b_head;
    while(b_temp->NEXT != NULL)
    {
        result += gaussian_chain_SIntegral(a_temp,b_temp);

        b_temp = b_temp->NEXT;
    }

    result += gaussian_chain_SIntegral(a_temp,b_temp);

    gaussian_chain_free(a_head);
    gaussian_chain_free(b_head);

    return result;
}

void orbital_S_matrix(gsl_matrix * dest, orbital * HEAD, int length)
{
    orbital * temp1, * temp2;

    temp1 = HEAD;
    temp2 = HEAD;

    int i,j;

    double temp;

    for(i=0;i<length;i++)
    {
        for(j=i;j<length;j++)
        {
            gsl_matrix_set(dest,i,j,orbital_SIntegral(orbital_enquiry(HEAD,i),orbital_enquiry(HEAD,j)));
        }
    }

    gsl_matrix_flip(dest,length,length);
}