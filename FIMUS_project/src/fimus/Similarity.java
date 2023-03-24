/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package fimus;
import java.util.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

/**
 *
 * @author Md Geaur Rahman
 *
 * implementing Helen's Similarity measure (S' (S-Prime) and S" (S-Secumdum))
 * Version 1.0 10/08/2012
 */
public class Similarity
{

    /*
     * Global declaration
     */
 private double [][]S_Prime;  //Similarity S_Prime
 private double [][]S_Secundum;  //Similarity S_Secundum
 private double [][]S;  //Weighted Similarity
 private double w1=0.6; //Weighted value for S_Prime //By default 0.6
 private double w2=1-w1; //Weighted value for S_Secundum
 private double Threshold=0.4; //Weighted value for calculating S_Secundum

 private int totalRecord;  //total records of the data file
 private int totalAttr;  //total records of the data file
 private String []domainValues;//contains domain values of all attributes
 private int []start;//strating position of each attribute
 private int []domainsize;//contains domain size of each attribute
 private int totalDomainValue;//contains total domain values of all attributes
 private int [][]C;  //CoAppearance matrix
 private String [][]dataset;//contains dataset
 private int []degree;//contains total degrees/edges of each category

 /**
     * calculate similarities among the categories of all attributes
     * @Param dataFile - contains data set
     * @return S - contains a similarity matrix
     */
 public double [][]similarityMeasure(String [][]dataFile)
  {
      dataset=dataFile;
      initialize(); //initilize all global variables
      CoAppearanceMatrix();// find coappearance matrix
      calculate_S_Prime();//calculate S-Prime
      calculate_S_Secundum();//calculate S-Secundum
      weightedSimilarity();//calculate weighted similarity
      normalizedSimilarity();// normalise
      return S;
  }
 /**
  * get co-apperance matrix
  */
 public int [][]getCAM()
 {
     return C;
 }
/**
  * get start location of each attribute
  */
 public int []getStart()
 {
     return start;
 }

 /**
  * get domain values
  */
 public String []getDominValues()
 {
     return domainValues;
 }
 /**
  * get total Domain Value
  */
 public int getTotalDomainValue()
 {
     return totalDomainValue;
 }
 /**
  * get DomainSize of each attribute
  */
 public int []getDomainSize()
 {
     return domainsize;
 }

/**
 * the method is used to initialize the global variables
 * @Param dataFile- contains data set
 */
 private void initialize()
   {
     
     totalRecord=dataset.length;
     totalAttr=dataset[0].length;
     start=new int[totalAttr];
     domainsize=new int[totalAttr];
     
      totalDomainValue=0;
     //find domain info of each attribute
     for(int i=0; i<totalAttr;i++)
     {
         String []tmpdomainValues=findDomainInfo(i);
         int tmpDS=totalDomainValue;
         String []tmpDV=new String[tmpDS];
         if(i!=0){
         System.arraycopy(domainValues, 0, tmpDV, 0, tmpDS);
         }
         start[i]=totalDomainValue;
         totalDomainValue+=domainsize[i];
         domainValues=new String[totalDomainValue];
         if(i==0)
         {
             System.arraycopy(tmpdomainValues, 0, domainValues, 0, totalDomainValue);
         }
         else
         {
             System.arraycopy(tmpDV, 0, domainValues, 0, tmpDS);
             System.arraycopy(tmpdomainValues, 0, domainValues, tmpDS, domainsize[i]);
         }
     }
     //initialize co-appearance matrix, similarity matrix
     S_Prime=new double[totalDomainValue][totalDomainValue];
     S_Secundum=new double[totalDomainValue][totalDomainValue];
     S=new double[totalDomainValue][totalDomainValue];
     C=new int[totalDomainValue][totalDomainValue];
     degree=new int[totalDomainValue];
     for(int i=0; i<totalDomainValue;i++)
     {
         degree[i]=0;
         for(int j=0; j<totalDomainValue;j++)
         {
             S_Prime[i][j]=0.0;S_Secundum[i][j]=0.0;S[i][j]=0.0;
             C[i][j]=0;
         }
     }
    }

 /**
     * calculate similarities among the categories of all attributes
     * @Param dataFile - contains data set
     * @return S - contains a similarity matrix
     */
 public double [][]similarityMeasure(String [][]dataFile,String [][]data, 
         int[]attrNtype,int[][]MV,int[]MR, double isSim)
  {
      dataset=dataFile;
      int tr=data.length;int tm=0;
      int attr=attrNtype.length;
      for(int i=0;i<tr;i++)
      for(int j=0;j<attr;j++)
      if(MV[i][j]==1)tm++;
      Format fm=new Format();
      if(isSim>=0.8&& isSim<=0.9)
        {
            fm.runFormat(data, MV, MR, attrNtype, tm);
            Initialize(); //initilize all global variables
        }
      else{
             initialize(); //initilize all global variables
        }
      CoAppearanceMatrix();// find coappearance matrix
      calculate_S_Prime();//calculate S-Prime
      calculate_S_Secundum();//calculate S-Secundum
      weightedSimilarity();//calculate weighted similarity
      normalizedSimilarity();// normalise
      return S;
  }

 /**
 * the method is used to initialize the global variables
 * @Param dataFile- contains data set
 */
 private void Initialize()
   {

     totalRecord=dataset.length;
     totalAttr=dataset[0].length;
     start=new int[totalAttr];
     domainsize=new int[totalAttr];

      totalDomainValue=0;
     //find domain info of each attribute
     for(int i=0; i<totalAttr;i++)
     {
         String []tmpdomainValues=FindDomainInfo(i);
         int tmpDS=totalDomainValue;
         String []tmpDV=new String[tmpDS];
         if(i!=0){
         System.arraycopy(domainValues, 0, tmpDV, 0, tmpDS);
         }
         start[i]=totalDomainValue;
         totalDomainValue+=domainsize[i];
         domainValues=new String[totalDomainValue];
         if(i==0)
         {
             System.arraycopy(tmpdomainValues, 0, domainValues, 0, totalDomainValue);
         }
         else
         {
             System.arraycopy(tmpDV, 0, domainValues, 0, tmpDS);
             System.arraycopy(tmpdomainValues, 0, domainValues, tmpDS, domainsize[i]);
         }
     }
     //initialize co-appearance matrix, similarity matrix
     S_Prime=new double[totalDomainValue][totalDomainValue];
     S_Secundum=new double[totalDomainValue][totalDomainValue];
     S=new double[totalDomainValue][totalDomainValue];
     C=new int[totalDomainValue][totalDomainValue];
     degree=new int[totalDomainValue];
     for(int i=0; i<totalDomainValue;i++)
     {
         degree[i]=0;
         for(int j=0; j<totalDomainValue;j++)
         {
             S_Prime[i][j]=0.0;S_Secundum[i][j]=0.0;S[i][j]=0.0;
             C[i][j]=0;
         }
     }
    }

 /**
 * calculate S-Prime
 * S-Prime_pq=(sum^N_k=1(sqrt(a_pk * a_kq))/sqrt(d_p * d_q)
 */
private void calculate_S_Prime()
    {
     
     for(int i=0; i<totalAttr;i++)
     {
         for(int p=start[i]; p<(start[i]+domainsize[i]);p++)
         {
            for(int q=start[i]; q<(start[i]+domainsize[i]);q++)
            {
                if(p==q)
                {
                    S_Prime[p][q]=1.0;
                }
                else
                {
                     double numerator=0.0,denumetor=0.0;
                     for(int k=0; k<totalDomainValue;k++)
                     {
                        if(k<=start[i]|| k>=(start[i]+domainsize[i]))
                        {
                            numerator+=Math.sqrt(C[p][k]*C[k][q]);
                        }
                     }
                     denumetor=Math.sqrt(degree[p]*degree[q]);

                     if(denumetor>0.0) S_Prime[p][q]=numerator/denumetor;
                }
             }
         }
     }
 }
/**
 * calculate S-Secundum
 * S-Secundum_pq=S-Prime_pq of the neighbours
 */
private void calculate_S_Secundum()
    {
    /* Loop over all attributes in the data set*/
     for(int x=0; x<totalAttr;x++)
     {
         /* Loop over all pairs of values for the attribute */
         for(int i=start[x]; i<(start[x]+domainsize[x]);i++)
         {
            for(int j=start[x]; j<(start[x]+domainsize[x]);j++)
            {
                if(i==j)
                {
                    S_Secundum[i][j]=1.0;continue;
                }
             int currIndex=0;
             int []neighList=new int[totalDomainValue];
             for(int k=0; k<totalDomainValue;k++)
             {
                 neighList[k]=-1;
             }
             /* Loop over all attributes in the data set, excluding x */
             for(int y=0; y<totalAttr;y++)
             {
                 if(x!=y)
                 {
                     /* Loop over all pairs of values in y */
                     for(int c=start[y]; c<(start[y]+domainsize[y]);c++)
                     {
                         for(int d=start[y]; d<(start[y]+domainsize[y]);d++)
                         {
                             if(S_Prime[c][d]>Threshold &&
                                     ((C[c][i]>0 &&C[d][j]>0)||(C[d][i]>0 &&C[c][j]>0)))
                             {
                                 //implementing function mergeNeighbours
                                 if(neighList[c]==-1 && neighList[d]==-1 )
                                 {
                                    neighList[c]=currIndex;
                                    neighList[d]=currIndex;
                                    currIndex++;
                                 }
                                 else if(neighList[c] != -1 && neighList[d] != -1)
                                 {
                                    if(neighList[c]!=neighList[d])
                                    {
                                        neighList[d]=neighList[c];
                                        for(int z=0; z<totalDomainValue;z++)
                                        {
                                            if(neighList[z]==neighList[d])
                                            {
                                                neighList[z]=neighList[c];
                                            }
                                        }
                                    }
                                 }
                                 else if(neighList[c] != -1 && neighList[d] == -1)
                                 {
                                     neighList[d]=neighList[c];
                                 }
                                 else if(neighList[c] == -1 && neighList[d] != -1)
                                 {
                                     neighList[c]=neighList[d];
                                  }

                             }
                         }

                     }

                 }
             }
             //calculating S-Prime on the merged neighbours (categories).
             //After merging, we have "currIndex" number of categories
             double numerator=0.0,denumetor=0.0;
             for(int z=0; z<currIndex;z++)
             {
                int tmpiT=0, tmpjT=0;
                for(int k=0; k<totalDomainValue;k++)
                {
                    if(neighList[k]==z)
                    {
                        tmpiT+=C[i][k];
                        tmpjT+=C[j][k];
                    }
                }
                numerator+=Math.sqrt(tmpiT*tmpjT);
             }
             denumetor=Math.sqrt(degree[i]*degree[j]);
             //calculating S_Secundum[i][j] for i-th and j-th categories
             if(denumetor>0.0) S_Secundum[i][j]=numerator/denumetor;
            }
         }
     }
 }

/**
 * calculate Weighted similarity
 * S_pq=S-Prime_pq * w1+ S-Secundum *w2
 */
private void weightedSimilarity()
{
     for(int p=0; p<totalDomainValue;p++)
     {
         for(int q=0; q<totalDomainValue;q++)
             {
                S[p][q]=S_Prime[p][q]*w1+S_Secundum[p][q]*w2;
             }
      }
 }
/*
 * Calculate normalized similarity
 */
private void normalizedSimilarity()
{
     for(int p=0; p<totalDomainValue;p++)
     {
         double total=0.0;
         for(int q=0; q<totalDomainValue;q++)
             {
                total+=S[p][q];
             }
         if(total>0)
         {
         for(int q=0; q<totalDomainValue;q++)
             {
                S[p][q]=S[p][q]/total;
             }
         }
      }
 }

/**
 * Find co-appearance matrix and degrees of a data set
 */
private void CoAppearanceMatrix()
    {
     //find coappearnce matrix for the data set
     for(int i=0; i<totalAttr-1;i++)
     {
         for(int p=start[i]; p<(start[i]+domainsize[i]);p++)
         {
            for(int j=i+1; j<totalAttr;j++)
            {
                for(int q=start[j]; q<(start[j]+domainsize[j]);q++)
                {
                    C[p][q]=findCoappearance(i,domainValues[p],j,domainValues[q]);
                    C[q][p]=C[p][q];
                }
            }
         }
     }
     //find degrees
     for(int i=0; i<totalDomainValue;i++)
     {
         for(int j=0; j<totalDomainValue;j++)
            {
             degree[i]+=C[i][j];
            }
      }

 }
/**
 * The method is used to find the co-appearance of x and y
 * @param
 * attr1-contains attribute index for x.
 * x- is a value of attr1
 * attr2-contains attribute index for y.
 * y- is a value of attr2
 * @return f- is the co-appearance of x and y.
 */
private int findCoappearance(int attr1,String x, int attr2,String y)
{
    int f=0;
    for(int i=0; i<totalRecord;i++)
    {
       if(dataset[i][attr1].equals(x) && dataset[i][attr2].equals(y))
       {
            f++;
        }
    }
    return f;
}
/**
 * The method is used to find domain information of a given attribute
 * @param
 * attrPos-contains attribute index for which domain info will be identified.
 * @return tmpDomain- contains domain values of the attribute.
 */
 private String []findDomainInfo(int attrPos)
    {
        int attrDomainSize=0;
        String []tmpDomain=new String[totalRecord];
        for(int i=0; i<totalRecord;i++)
        {
            if(isMissing(dataset[i][attrPos])==0)
            {
               if(chkDomain(tmpDomain,attrDomainSize,dataset[i][attrPos])==0)
                {
                    tmpDomain[attrDomainSize]=dataset[i][attrPos];
                    attrDomainSize++;
                }
            }
        }
        domainsize[attrPos]=attrDomainSize;
        return tmpDomain;

    }

 /**
 * The method is used to find domain information of a given attribute
 * @param
 * attrPos-contains attribute index for which domain info will be identified.
 * @return tmpDomain- contains domain values of the attribute.
 */
 private String []FindDomainInfo(int attrPos)
    {
        int attrDomainSize=0;
        String []tmpDomain=new String[totalRecord];
        for(int i=0; i<totalRecord;i++)
        {
            if(isMissing(dataset[i][attrPos])==0)
            {
               if(chkDomain(tmpDomain,attrDomainSize,dataset[i][attrPos])==0)
                {
                    tmpDomain[attrDomainSize]=dataset[i][attrPos];
                    attrDomainSize++;
                }
            }
        }
        domainsize[attrPos]=0;
        return tmpDomain;

    }
 /**
 * The method is used to check whether or not a given value is already in the
 * domain list.
 * @param
 * tmpDomain-contains domain values of an attribute.
 * domainSize-the total number of values of the attribute
 * curVal- the current value is to be checked with existing domain values.
 * @return flag- is an integer value indicating Exist (1) or NOT exist(0)
 */
 private int chkDomain(String []tmpDomain,int domainSize, String curVal)
    {
        int flag=0;
        for(int i=0;i<domainSize;i++)
        {
            if(curVal.equals(tmpDomain[i]))
            {
               flag=1; break;
            }
        }
        return flag;
    }
 /**
  * the method prints the values of an given array
  * @param
  * Msg- the message to be printed as a label.
  * sv-contains similarity values of the attribute categories
  */
private void printArray(String Msg, double [][]sv)
{
        System.out.print(Msg+"\n");
        DecimalFormat df = new DecimalFormat();
        df.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
        df.setMinimumIntegerDigits(1);
        df.setMaximumFractionDigits(3);
        df.setMinimumFractionDigits(3);
        df.setGroupingUsed(false);
        for(int i=0;i<sv.length;i++)
        {
            for(int j=0;j<sv[0].length;j++)
            {
                System.out.print(df.format(sv[i][j])+", ");
            }
            System.out.print("\n");
         }
}
/**
  * the method prints the values of an given array
  * @param
  * Msg- the message to be printed as a label.
  * sv-contains similarity values of the attribute categories
  */
private void printIntArray(String Msg, int [][]sv)
{
        System.out.print(Msg+"\n");
        for(int i=0;i<sv.length;i++)
        {
            for(int j=0;j<sv[0].length;j++)
            {
                System.out.print(sv[i][j]+",  ");
            }
            System.out.print("\n");
         }
}
/**
  * this function will indicate whether or not a value is missing.
  *
  * @param oStr the string to be checked
  * @return ret an integer value 0->No missing, 1->Missing
  */

 private int isMissing(String oStr)
    {
       int ret=0;
       if(oStr.equals("")||oStr.equals("?")||oStr.equals("�"))
                     {
                         ret=1;
                    }
       return ret;
    }
 


}
