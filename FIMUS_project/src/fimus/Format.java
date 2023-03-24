/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package fimus;

import java.io.*;
import Jama.*;
import java.util.Random;
/**
 * A logic rule {@link Node} of a {@link DecisionTree}. This node has the additional
 * information relating to class supports. That is, for the data set partition
 * of records that 'fall' into this leaf node, count the number of records having
 * each class value. Also has information about the class itself.
 *
 * @author helengiggins
 * @version 1.0 30/11/2010
 * @see Node
 * @see DecisionTree
 */
public class Format
{
    /** used for managing reading and writing to file */

double[] attrMean;
String  []aty;
/**

 
/**
 * 
 *
 * @param data 2D array 
 * @param FormatMatrix 2D array
 * @param FormatRecord 1D array
 * @param attrNType 1D array where 1-> numerical, 0->categorical
 */
    
public void runFormat(String [][]data, int [][]FormatMatrix,
        int[] FormatRecord, int[] attrNType, int totalRule)
{
    int totalRecords=data.length;
    int totalAttrs=data[0].length;
    attrMean=new double[totalAttrs];


    int numAttr=0;
    for(int c=0; c<totalAttrs;c++)
    {
            if(attrNType[c]==1)
            {
                numAttr++;
            }
    }
   
    
//    if(totalRecords>numAttr)
//    {
        formatalgorithm(data,totalRecords, totalAttrs,numAttr,
             totalRule, FormatMatrix, FormatRecord, attrNType, 1);
//    }


}


private void calAttrMean(String [][]data, int []attrType)
{
    int totalAttrs=data[0].length;
    int totalRecords=data.length;
    for(int c=0; c<totalAttrs;c++)
        {
            if(attrType[c]==1)
            {
                double total=0.0; int cnt=0;
                for(int i=0; i<totalRecords;i++)
                {
                    if(isMissing(data[i][c])==0)
                    {
                        total+=Double.parseDouble(data[i][c]);
                        cnt++;
                    }
                }
                if(cnt>0)attrMean[c]=total/(double)cnt;
            }

        }
}

/**
  * this function 
  *
  * @param oStr the string to be checked
  * @return ret 
  */

 private int isMissing(String oStr)
    {
       int ret=0;
       if(oStr.equals("")||oStr.equals("?")||oStr.equals("�")||oStr.equals("NaN")||oStr.equals("  NaN"))
                     {
                         ret=1;
                    }
       return ret;
    }

 

    /**
  * set attr info and return no. of numerical attr
  */
    private int getAttrType(String attrFile,int []attrNType,String  []attrSType)
    {
         FileManager fileManager=new FileManager();
         String [][]tmpAty=fileManager.readFileAs2DArray(new File(attrFile));
         int nAttr=tmpAty[0].length;
         int numAttr=0;
         for(int i=0; i<nAttr;i++)
         {
             if(tmpAty[0][i].equals("1"))
             {
                 attrNType[i]=1;
                 attrSType[i]="n";
                 aty[i]="n";
                 numAttr++;
             }
            else
             {
                 attrNType[i]=0;
                 attrSType[i]="c";
                 aty[i]="c";
             }
         }
        return numAttr;
    }
    
    private int chkAvg(String [][]dFile,int noOfRecords,int noOfAttrs,
            int []attrType,int sPos,int []avgImputation,
            int []NewAttrNType,double[] attrAvgVal )
    {
         
         int totAvg=0;
         double preVal, curVal;
         String str="";
         for(int c=sPos;c<noOfAttrs;c++)
         {
             if (attrType[c]==1)
            {
                int cnt=0;
                preVal=0.0;curVal=0.0;
                for(int i=0; i<noOfRecords;i++)
                {
                    str=dFile[i][c];
                    if (isMissing(str)==1)
                        curVal=0.0;
                    else
                        curVal=Double.parseDouble(str);
                    if(curVal!=preVal)
                    {
                        preVal=curVal;
                        cnt++;
                    }

                   if(cnt>=2)break;
                }
                if(cnt<=1&& noOfRecords!=1)
                 {
                        avgImputation[c]=1;
                        NewAttrNType[c]=0;
                        attrAvgVal[c]=curVal;
                        totAvg++;
                 }
             }
         }
         return totAvg;
    }
    /*
     * Implementation of the format algorithm
     */
    public void formatalgorithm(String [][]data,int noOfRecords,
            int noOfAttrs,int noOfNumericalAttrs, int totalMissingValues,
            int [][]MissingMatrix, int[] MissingRecord,
            int []attrNType, int MissAvailable)
    {
        double [][]dataOriginal=new double[noOfRecords][noOfNumericalAttrs];

        double val;
        int [][]nMissingMatrix=new int[noOfRecords][noOfNumericalAttrs]; 
        int[] nMissingRecord=new int[noOfRecords];
        
        int nc=0;
        
        for(int j=0;j<noOfAttrs;j++)
        {
            if (attrNType[j]==1)
            {
                for(int i=0; i<noOfRecords;i++)
                {
                        if (MissingMatrix[i][j]==1)
                        {
                           if(isMissing(data[i][j])==1)
                               val = 0.0;
                           else
                               val=Double.parseDouble(data[i][j]);
                            nMissingMatrix[i][nc]=1;
                           
                        }
                        else
                        {
                           val=Double.parseDouble(data[i][j]);
                           nMissingMatrix[i][nc]=0;
                        }
                        dataOriginal[i][nc]=val;
                }
                nc++;
            }
        }
        for(int i=0; i<noOfRecords;i++)
        {
            nMissingRecord[i]=0;
            for(int j=0;j<noOfNumericalAttrs;j++)
            {
                if (nMissingMatrix[i][j]==1)
                 {
                    nMissingRecord[i]=1;break;
                }
            }

        }

        double [][]dataCurrent=new double[noOfRecords][noOfNumericalAttrs];
        double []mu=new double[noOfNumericalAttrs];
        initilizeArray(mu);
        double []muPrevious=new double[noOfNumericalAttrs];
        double [][]cov=new double[noOfNumericalAttrs][noOfNumericalAttrs];
        double [][]covPrevious=new double[noOfNumericalAttrs][noOfNumericalAttrs];
           
       MatrixCalculation mxCal=new MatrixCalculation();

       //copy of original data for updating
       for(int i=0; i<noOfRecords;i++)
       {
           for(int j=0; j<noOfNumericalAttrs;j++)
           {
                dataCurrent[i][j]=dataOriginal[i][j];
           }
        }
        int T=0;

        do{
           for(int i=0; i<noOfNumericalAttrs;i++)
           {
               muPrevious[i]=mu[i];
               for(int j=0; j<noOfNumericalAttrs;j++)
               {
                    covPrevious[i][j]=cov[i][j];
               }
            }
           
           computeMean(dataCurrent, nMissingMatrix, MissAvailable, mu);
           mxCal.computeCovariance(dataCurrent, noOfRecords, noOfNumericalAttrs,  mu,cov);

           for(int i=0;i<noOfRecords;i++)
           {
               if(nMissingRecord[i]==1)
               {
                Leaves(dataCurrent,i,nMissingMatrix,T,mu,cov,noOfNumericalAttrs);
               }
           }
           T++;
       }while(T<1);

       
       nc=0;
       for(int c=0;c<noOfAttrs;c++)
       {
           if(attrNType[c]==1)
           {
               for(int i=0; i<noOfRecords;i++)
               {
                   if(MissingMatrix[i][c]==1)
                   {
                       if(dataCurrent[i][nc]!=Double.NaN)
                       {
                            data[i][c]= dataCurrent[i][nc]+"";
                       }
                       if(isMissing(data[i][c])==1)
                       {
                            data[i][c]=attrMean[c]+"";
                       }
                   }
               }
               nc++;
           }
       }
    }
   

   
    public void Leaves(double [][]dataElement,int msRow,int [][]MissingMatrix,
            int T, double[]mu,double [][]cov,
            int noOfNumericalAttrs)
    {
       int r,c,k,j;
       int m=0;
       
       for(j=0; j<noOfNumericalAttrs;j++)
       {
           if(MissingMatrix[msRow][j]==1)m++;
        }

       int a=noOfNumericalAttrs-m;
       double []X_a=new double[a];
       double []X_m=new double[m];
       double []Mu_a=new double[a];
       double []Mu_m=new double[m];
       double [][]cov_aa=new double[a][a];
       double [][]cov_mm=new double[m][m];
       double [][]cov_am=new double[a][m];

      
       Matrix matX_a,matX_m, matMu_a,matMu_m,matcov_aa,matcov_mm,matcov_am,matcov_ma;
       Matrix matinvOfcov_aa,matB,matC,mate;
       k=0;
       for(c=0;c<noOfNumericalAttrs;c++)
       {
           if(MissingMatrix[msRow][c]==0)
           {
               X_a[k]=dataElement[msRow][c];
               Mu_a[k]=mu[c];
               k++;
           }
       }

       MatrixCalculation mxCal=new MatrixCalculation();
      
       matX_a= mxCal.oneDArrayToMatrix(X_a);
       matMu_a= mxCal.oneDArrayToMatrix(Mu_a);
        for(c=0,k=0;c<noOfNumericalAttrs;c++)
        {
            if(MissingMatrix[msRow][c]==1)
            {
                X_m[k]=dataElement[msRow][c];
                Mu_m[k]=mu[c];
                if (k<m)k++;
            }
        }
       
       matX_m= mxCal.oneDArrayToMatrix(X_m);
       matMu_m= mxCal.oneDArrayToMatrix(Mu_m);
       if(a>0)
       {
           k=0;
           for(r=0;r<noOfNumericalAttrs;r++)
           {
               j=0;
               if(MissingMatrix[msRow][r]==0)
               {
                   for(c=0;c<noOfNumericalAttrs;c++)
                   {
                       if(MissingMatrix[msRow][c]==0)
                       {
                           cov_aa[k][j]=cov[r][c];
                           j++;
                       }
                   }
                   k++;
               }
            }
           matcov_aa= new Matrix(cov_aa);
          double tt=matcov_aa.det();
          if(tt!=0)
          {
    
          k=0;
           for(r=0;r<noOfNumericalAttrs;r++)
           {
               j=0;
               if(MissingMatrix[msRow][r]==1)
               {
                   for(c=0;c<noOfNumericalAttrs;c++)
                   {
                       if(MissingMatrix[msRow][c]==1)
                       {
                           cov_mm[k][j]=cov[r][c];
                           j++;
                       }
                   }
                   k++;
               }
            }


           matcov_mm=new Matrix(cov_mm);
           k=0;
           for(r=0;r<noOfNumericalAttrs;r++)
           {
               j=0;
               if(MissingMatrix[msRow][r]==0)
               {
                   for(c=0;c<noOfNumericalAttrs;c++)
                   {
                       if(MissingMatrix[msRow][c]==1)
                       {
                           cov_am[k][j]=cov[r][c];
                           j++;
                       }
                   }
                   k++;
               }
            }
          matcov_am= new Matrix(cov_am);
          matcov_ma=matcov_am.transpose(); 
  
          try
          {
          matinvOfcov_aa=matcov_aa.inverse();
          matB=matinvOfcov_aa.times(matcov_am);
          matC=matcov_ma.times(matinvOfcov_aa.times(matcov_am));
          ////
          Matrix Xa_Minus_Mua=matX_a.minus(matMu_a);
          Matrix Xa_Minus_Mua_TimesB=Xa_Minus_Mua.times(matB);
          matX_m=matMu_m.plus(Xa_Minus_Mua_TimesB);
//          if(T==0)
//          {
//            mate = generateResidualMatrix(X_m,matMu_m, matcov_mm,matC, msRow, m);
//            matX_m=matX_m.plus(mate);
//            }
          }
          catch(Exception ex)
          {
            matX_m=matMu_m;
          }
           }
         else
          {
              matX_m=matMu_m;
            }
        }
       else
        {
           matX_m=matMu_m;
        }
      //result
      mxCal.matrixTo1DArray(matX_m,X_m,m);
      // update dataset
       for(c=0,k=0;c<noOfNumericalAttrs;c++)
        {
            if(MissingMatrix[msRow][c]==1)
            {
                dataElement[msRow][c]=X_m[k];
                if(k<m)k++;
            }
        }
      
    }
public Matrix  generateResidualMatrix(double []matX_m,Matrix matMu_m,
            Matrix matCov_m,Matrix pertVar,int recNum,int m)
    {
            Matrix output, cholPertVariance,pertVariance;
            Matrix mean, randomMatrix;
            CholeskyDecomposition cholesky;
            pertVariance=pertVar;
            cholesky = pertVariance.chol();
            cholPertVariance = cholesky.getL();

            mean = calculateMeanZero(matX_m,matMu_m,matCov_m,recNum,m);//return zero mean vector
            randomMatrix = generateRandomVariates(m);
            output = mean.plus(pertVar.times(randomMatrix.transpose()));
            output = output.transpose(); // it produces a matrix of dim.(1x6)
            return output;
    }

    public Matrix calculateMeanZero(double []X_m,Matrix matMu_m,Matrix matCov_m,int recNum,int m){
        Matrix c, mean, result,covXX,covYX,covYY;
        MatrixCalculation mxCal=new MatrixCalculation();

        double [] c_rec = new double[m];
        for (int i=0;i<m;i++)
            c_rec[i]=0.0;//X_m[i];
            c = mxCal.oneDArrayToMatrix(c_rec);
        return c.transpose();//mean;
    }

    public Matrix generateRandomVariates(int numAttr){
        Matrix rndMat;
        double[] randVariates = new double[numAttr];
        for(int i=0; i<numAttr; i++)
           randVariates[i] = 0;
        Random rand = new Random();
        for(int j=0; j<numAttr; j++){
           randVariates[j] = rand.nextGaussian();
        }
        MatrixCalculation mxCal=new MatrixCalculation();
        rndMat =  mxCal.oneDArrayToMatrix(randVariates);
       return rndMat;
    }


    /*
     * Calculate mean based on the available values.
     */
    private void computeMean(double[][]DataElement, int [][]MissingMatrix,int msFlag,double []mu )
    {
        int totRec;
        int numAttr=DataElement[0].length;
        int numRecords=DataElement.length;
        for(int i=0; i<numAttr; i++)
            mu[i] = 0;
        double sum=0;
        for(int i=0; i<numAttr; i++)
        {
           sum = 0.0;totRec=0;
           for(int j = 0; j<numRecords; j++)
           {
                   
                   if(MissingMatrix[j][i]==0 ||msFlag==1) 
                   {
                       totRec++;
                       sum += DataElement[j][i];
                   }
           }
          if(totRec>0)     mu[i] = sum/(double)totRec;
        }

    }

    //check whether an attribute has found in the array msCol.
public int chkAttrIsMiss(int sCol,int []msCol,int sPos,int ePos)
    {
        int i, flag=0;
        for(i=sPos;i<ePos;i++)
        {
            if(msCol[i]==sCol)
            {
                flag=1;break;
            }
        }
        return flag;
    }

   
    public int calMeanError(double []preMu,double []curMu)
    {
       
        int isError=1, tA=0;
        double diff=0.0;
        double cuT=0.0,cyAvg=0.0;
        double preT=0.0,preAvg=0.0;
        tA=preMu.length;
        for(int i=0;i<tA;i++)
        {
            cuT=cuT+curMu[i];
            preT=preT+preMu[i];
        }
        if(tA>0)
        {
            cyAvg= cuT/tA;
            preAvg= preT/tA;
        }

        diff=Math.abs(preAvg-cyAvg);

         if(diff<0.0000000001)
         {
            isError=0;
         }

        
        return isError;
    }

  
    public int calCovError(double [][]preCov,double [][]curCov)
    {

        int isError=1;
        double diff=0.0;
        double preDet=0.0;
        double curDet=0.0;
        try
        {
            Matrix preCov_mat=new Matrix(preCov);
            Matrix curCov_mat=new Matrix(curCov);
            preDet=preCov_mat.det();
            curDet=curCov_mat.det();
            diff=Math.abs(preDet-curDet);
            if(diff<0.0000000001)
             {
                isError=0;
             }
        }
        catch(Exception e)
        {
         }
        return isError;
    }

    
    //following method will initialize an array.(data type double)  
     public void initilizeArray(double []curArr)
    {
      for(int i=0;i<curArr.length;i++)
          curArr[i]=0.0;
     }

     
}
