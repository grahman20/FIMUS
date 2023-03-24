/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package fimus;
import java.io.*;
import java.util.*;
import java.text.DecimalFormat;

/**
 * FIMUS imputes numerical and categorical missing values by using a data set’s 
 * existing patterns including co-appearances of attribute values, correlations 
 * among the attributes and similarity of values belonging to an attribute. 
 * 
 * <h2>Reference</h2>
 * 
 * Rahman, M. G. and Islam, M. Z. (2014): FIMUS: A Framework for Imputing Missing Values Using Co-Appearance, Correlation and Similarity Analysis, Knowledge-Based Systems, 56, 311 - 327, ISSN 0950-7051, DOI information: http://dx.doi.org/10.1016/j.knosys.2013.12.005
 *  
 * @author Md Geaur Rahman <https://csusap.csu.edu.au/~grahman/>
 */

/**
 *
 * @author grahman
 * 10/8/2012
 */
public class FIMUS
{

    /*
     * Global declaration
     */
 public double uThreshold; //weight on the values of the record contains missing values
 private double lamda; //weight on the attribute similar values of the record contains missing values

 private int totalRecord;  //total records of the data file
 private int totalAttr;  //total records of the data file
 private String []domainValues;//contains domain values of all attributes
 private int []start;//strating position of each attribute
 private int []domainsize;//contains domain size of each attribute
 private int totalDomainValue;//contains total domain values of all attributes
 private int [][]C;  //CoAppearance matrix
 private int [][]MV;  //Maissing values, 0->no missing, 1->Missing
 private int []MR;  //Missing Records, 0->no missing, 1->Missing
 private long totalMissing;
 private double [][]R;  //Correlation matrix
 private int [][]ATotal;  //attribute total
 private double iSim;//contains formatiing information;
 private int [][]RowTotal;  //Row total
 private int [][]ColTotal;  //Column total
 private int t=0;
 private String [][]dataset;//contains dataset
 private String [][]datasetG;//contains generalized dataset
 private double [][]S;//contains total degrees/edges of each category
 private int [] attributeType; //contain attributes type 0->Categorical, 1->Numerical
 private String [] attrSType; //contain attributes type c->Categorical, n->Numerical
     /**
     * Impute a given data set having missing values
     * and calls other necessary methods
     *
     * @param attrFile contains 2 lines attributes types and name information
     * @param dataFile data set having missing values to be imputed
     * @param outputFile filename of the imputed data set
     */
 public void runFIMUS(String attrFile, String dataFile,String outputFile, double userThreshold)
 {
    uThreshold=userThreshold;
    lamda=1-uThreshold;

    FileManager fileManager = new FileManager();

    dataset=fileManager.readFileAs2DArray(new File(dataFile));
    totalRecord=dataset.length;
    totalAttr=dataset[0].length;
    datasetG=new String[totalRecord][totalAttr];
    copy2DArray(dataset,datasetG);
    findMissing();
    getAttrType(attrFile);
    generalise();
    initialize(datasetG); //initialize the global variables


    cmiImpute(attrFile,dataFile);
    iterative(attrFile,dataFile);
    arrayToFile(dataset,totalRecord,totalAttr,outputFile);
 }
 //iterative imputation
 private void iterative(String attrF, String dataF)
 {
    int flag=1;
    double rmse=1.0, rmsep=1.0;
    do
    {
        String [][]dataPre=new String[totalRecord][totalAttr];
        copy2DArray(dataset,dataPre);
        copy2DArray(dataset,datasetG);

        generalise();
        initialize(datasetG); //initialize the global variables
        cmiImpute(attrF,dataF);
        rmse=calculateRMSE(dataPre,dataset,MV, MR, attributeType);
        if(rmse==0.0 || (rmse==rmsep))
            flag=0;
        rmsep=rmse;
    }while(flag==1);
 }
 //impute all missing values
  private void cmiImpute(String attrF, String dataF)
  {

//     FileManager fileManager = new FileManager();
//     String [][]tmattrF=fileManager.readFileAs2DArray(new File(attrF));
//     String [][]tmdataset=fileManager.readFileAs2DArray(new File(dataF));
     Similarity sm=new Similarity();

    //impute all values by cosiering as categorical values
    for(int i=0; i<totalRecord;i++)
     {
       if(MR[i]==1)
       {
        for (int j = 0; j < totalAttr; j++)
         {
            if(MV[i][j]==1)
            {
                datasetG[i][j]=CSR( datasetG[i], j, totalAttr, start, domainsize,domainValues,
          RowTotal, C, S,R);
               if( attributeType[j]==0)dataset[i][j]=datasetG[i][j];
            }
        }
      }
    }
    //impute all numerical missing  values
    for(int i=0; i<totalRecord;i++)
     {
       if(MR[i]==1)
       {
        for (int j = 0; j < totalAttr; j++)
         {
            if(MV[i][j]==1 && attributeType[j]==1)
            {
                String val=datasetG[i][j];
                int totalRec=0;
                int []tmp=new int[totalRecord];
                for(int k=0; k<totalRecord;k++)
                {
                    if(val.equals(datasetG[k][j]))
                    {
                        tmp[totalRec]=k;totalRec++;
                    }
                }
               String [][]datasetN=new String[totalRec][totalAttr];
               for(int k=0; k<totalRec;k++)
                {
                  System.arraycopy(dataset[tmp[k]], 0, datasetN[k], 0, totalAttr);
                }

                 double [][]S1=sm.similarityMeasure(datasetN,
                         dataset,attributeType,MV,MR,iSim);
                 int []start1=sm.getStart();
                 int []domainsize1=sm.getDomainSize();
                 String []domainValues1=sm.getDominValues();
                 if(domainsize1[j]>1)//calculate voting for more domain values
                 {
                     int totalDomainValue1=sm.getTotalDomainValue();
                     int [][]C1=sm.getCAM();


                     int [][]ATotal1=new int[totalAttr][totalAttr];
                     int [][]RowTotal1=new int[totalAttr][totalDomainValue1];
                     int [][]ColTotal1=new int[totalAttr][totalDomainValue1];
                     calRCTotal(totalAttr,start1,domainsize1,C1,ATotal1,RowTotal1,ColTotal1);

                     double [][]R1=new double[totalAttr][totalAttr];
                     calCorrelationV(totalAttr,start1,domainsize1,C1,ATotal1,RowTotal1,ColTotal1,totalDomainValue1, R1);
                     dataset[i][j]=CSR( dataset[i], j, totalAttr, start1, domainsize1,domainValues1,
              RowTotal1, C1, S1,R1);
                }
                 else if(domainsize1[j] == 1)
                 {
                        dataset[i][j]=domainValues1[start1[j]];
                 }
               }
        }
      }
//     arrayToFile(dataset,totalRecord,totalAttr,dataF);

    }
  }
 /*
  * impute a missing value
  */
 private String CSR( String []record, int ma, int totatr, int []srt, int []ds,String []dv,
          int [][]Rtotal, int [][]CAM, double [][]Sim,double [][]Corel)
 {
     int k=0;
     double []cv=new double [ds[ma]];
     for(int x=srt[ma]; x<srt[ma]+ds[ma];x++,k++)
     {
         cv[k]=0.0;
         for(int p=0; p<totatr;p++)
         {
             if (p!=ma && isMissing(record[p])==0 && Corel[ma][p]>0)
             {
                 double weightedVote=0.0;
                 double weightedVoteWithSimilarity=0.0;
                 double weightedVoteWithoutSimilarity=0.0;
                 int z=Integer.MAX_VALUE;
                 for(int y=srt[p]; y<srt[p]+ds[p];y++)
                 {
                     if(dv[y].equals(record[p]))
                     {
                         z=y;break;
                     }
                 }

                 if(uThreshold>0 && Rtotal[p][x]>0)
                 {
                    weightedVoteWithoutSimilarity = (double)CAM[x][z] / (double)Rtotal[p][x];
                 }

                 if(uThreshold<1.0)
                 {
                     for(int y=srt[p]; y<srt[p]+ds[p];y++)
                     {
                         double vote=0.0;
                         if(Rtotal[p][x]>0)
                         {
                             vote = (double)CAM[x][y] / (double)Rtotal[p][x];
                         }

                         weightedVoteWithSimilarity+=vote*Sim[z][y];
                     }
                 }
                 weightedVote=weightedVoteWithoutSimilarity*uThreshold+weightedVoteWithSimilarity*lamda;

                 cv[k]+=weightedVote*Corel[ma][p];
             }
        }

     }
    double cvM=Double.NEGATIVE_INFINITY;int MI=Integer.MAX_VALUE;
    for(k=0; k<ds[ma];k++)
    {
        if(cv[k]>cvM)
        {
            cvM=cv[k];
            MI=k;
        }
    }

    MI+=srt[ma];
    return dv[MI];
 }

 /*
  * impute a missing value
  */
 private String CSRi( String []record, int ma, int totatr, int []srt, int []ds,String []dv,
          int [][]Rtotal, int [][]CAM, double [][]Sim,double [][]Corel)
 {
     int k=0;
     double []cv=new double [ds[ma]];
     for(int x=srt[ma]; x<srt[ma]+ds[ma];x++,k++)
     {
         cv[k]=0.0;
         for(int p=0; p<totatr;p++)
         {
             if (p!=ma && isMissing(record[p])==0 && Corel[ma][p]>0)
             {
                 double weightedVote=0.0;
                 double weightedVoteWithSimilarity=0.0;
                 double weightedVoteWithoutSimilarity=0.0;
                 int z=Integer.MAX_VALUE;
                 for(int y=srt[p]; y<srt[p]+ds[p];y++)
                 {
                     if(dv[y].equals(record[p]))
                     {
                         z=y;break;
                     }
                 }

                 if(uThreshold>0 && Rtotal[p][x]>0)
                 {
                    weightedVoteWithoutSimilarity = (double)CAM[x][z] / (double)Rtotal[p][x];
                 }

                 if(uThreshold<1.0)
                 {
                     for(int y=srt[p]; y<srt[p]+ds[p];y++)
                     {
                         double vote=0.0;
                         if(Rtotal[p][x]>0)
                         {
                             vote = (double)CAM[x][y] / (double)Rtotal[p][x];
                         }

                         weightedVoteWithSimilarity+=vote*Sim[z][y];
                     }
                 }
                 weightedVote=weightedVoteWithoutSimilarity*uThreshold+weightedVoteWithSimilarity*lamda;

                 cv[k]+=weightedVote*Corel[ma][p];
             }
        }

     }
    double cvM=Double.NEGATIVE_INFINITY;int MI=Integer.MAX_VALUE;
    for(k=0; k<ds[ma];k++)
    {
        if(cv[k]>cvM)
        {
            cvM=cv[k];
            MI=k;
        }
    }

    MI+=srt[ma];
    return dv[MI];
 }


 /**
  * return attr info
  */
    private void getAttrType(String attrFile)
    {
     FileManager fileManager=new FileManager();
     String [][]tmpAty=fileManager.readFileAs2DArray(new File(attrFile));
     attributeType=new int[totalAttr];
     attrSType=new String[totalAttr];
     for(int i=0; i<totalAttr;i++)
     {
         if(tmpAty[0][i].equals("1"))
         {
             attributeType[i]=1;
             attrSType[i]="n";
         }
        else
         {
             attributeType[i]=0;
             attrSType[i]="c";
         }
     }

    }
 //find missing values
  private void findMissing()
  {
    int flg;
    totalMissing=0;
    MV=new int[totalRecord][totalAttr];
    MR=new int[totalRecord];
    for(int i=0; i<totalRecord;i++)
     {
        flg=0;
        for (int j = 0; j < totalAttr; j++)
         {
            MV[i][j]=isMissing(dataset[i][j]);
            if(MV[i][j]==1) {flg=1;totalMissing++;}
        }
        MR[i]=flg;
    }
  }
 /**
 * the method is used to initialize the global variables
 * @Param dataFile- contains data set
 */
 private void initialize(String [][]dataFile)
   {
     Similarity sm=new Similarity();
     S=sm.similarityMeasure(dataFile);
     start=sm.getStart();
     domainsize=sm.getDomainSize();
     totalDomainValue=sm.getTotalDomainValue();
     C=sm.getCAM();
     domainValues=sm.getDominValues();

     ATotal=new int[totalAttr][totalAttr];
     RowTotal=new int[totalAttr][totalDomainValue];
     ColTotal=new int[totalAttr][totalDomainValue];
     calRCTotal(totalAttr,start,domainsize,C,ATotal,RowTotal,ColTotal);
     R=new double[totalAttr][totalAttr];iSim=lamda;
     
     calCorrelationV(totalAttr,start,domainsize,C,ATotal,RowTotal,ColTotal,totalDomainValue, R);
    }

 //calculate row total and column total
 private void calRCTotal(int totatr, int []srt, int []ds,int [][]CAM,
         int [][]Atotal, int [][]Rtotal, int [][]Ctotal)
 {
     for(int i=0; i<totatr;i++)
     {
         for (int j = 0; j < totatr; j++)
         {
             //calculate row total
             int ijTot=0;
              for(int x=srt[j]; x<srt[j]+ds[j];x++)
              {
                int tmp=0;
                 for(int y=srt[i]; y<srt[i]+ds[i];y++)
                 {
                     tmp+=CAM[x][y];
                 }
                 Rtotal[i][x]=tmp;
                 ijTot+=tmp;
             }
             Atotal[i][j]=ijTot;

              //calculate column total
             for(int y=srt[j]; y<srt[j]+ds[j];y++)

              {
                int tmp=0;
                 for(int x=srt[i]; x<srt[i]+ds[i];x++)
                 {
                     tmp+=CAM[x][y];
                 }
                 Ctotal[i][y]=tmp;
             }
         }
     }
 }
 /*
  * Find correlation using Pearson's contingency co-efficient
  */
 private void calCorrelation(int totatr, int []srt, int []ds,int [][]CAM,
         int [][]Atotal, int [][]Rtotal, int [][]Ctotal,int tdv, double [][]Corel)
 {
     double [][]exCAM=new double [tdv][tdv];
     for(int i=0; i<totatr-1;i++)
     {
         for (int j = i+1; j < totatr; j++)
         {

             if (Atotal[i][j]>0)
             {
                 //calculate expected CAM
                  for(int x=srt[i]; x<srt[i]+ds[i];x++)
                  {
                     for(int y=srt[j]; y<srt[j]+ds[j];y++)
                     {
                         exCAM[x][y]=((double)(Rtotal[j][x]*Ctotal[i][y]))/((double)(Atotal[i][j]));
                     }
                  }
                  //calculate chi-sqaure
                  double chi=0.0;
                  for(int x=srt[i]; x<srt[i]+ds[i];x++)
                  {
                     for(int y=srt[j]; y<srt[j]+ds[j];y++)
                     {
                        if(exCAM[x][y]>0)
                            chi+=Math.pow(CAM[x][y]-exCAM[x][y],2)/exCAM[x][y];
                     }
                  }
                 Corel[i][j] = Math.sqrt(chi/(chi+Atotal[i][j]));
                 Corel[j][i]=Corel[i][j];
             }
            else
             {
                 Corel[i][j] = 0.0;
            }

         }
     }
 }
/*
  * Find correlation using Cramer's co-efficient
  */
 private void calCorrelationV(int totatr, int []srt, int []ds,int [][]CAM,
         int [][]Atotal, int [][]Rtotal, int [][]Ctotal,int tdv, double [][]Corel)
 {
     double [][]exCAM=new double [tdv][tdv];
     for(int i=0; i<totatr-1;i++)
     {
         for (int j = i+1; j < totatr; j++)
         {

             if (Atotal[i][j]>0)
             {
                 //calculate expected CAM
                  for(int x=srt[i]; x<srt[i]+ds[i];x++)
                  {
                     for(int y=srt[j]; y<srt[j]+ds[j];y++)
                     {
                         exCAM[x][y]=((double)(Rtotal[j][x]*Ctotal[i][y]))/((double)(Atotal[i][j]));
                     }
                  }
                  //calculate chi-sqaure
                  double chi=0.0;
                  for(int x=srt[i]; x<srt[i]+ds[i];x++)
                  {
                     for(int y=srt[j]; y<srt[j]+ds[j];y++)
                     {
                        if(exCAM[x][y]>0)
                            chi+=Math.pow(CAM[x][y]-exCAM[x][y],2)/exCAM[x][y];
                     }
                  }
                 int k=ds[i];
                 if(ds[j]<k)k =ds[j];
                 if (k<2) k=2;
                 Corel[i][j] = Math.sqrt(chi/(Atotal[i][j]*(k-1)));
                 Corel[j][i]=Corel[i][j];
             }
            else
             {
                 Corel[i][j] = 0.0;
            }

         }
     }
 }


 /*
 * this method is used to generalise all numerical attributes
 * into sqrt|domainsize| categories
 */
private void generalise()
    {
        int i,j;

        for(int k=0; k<totalAttr;k++)
        {
         if(attributeType[k]==1)
         {
            double max=Double.NEGATIVE_INFINITY;
            double min=Double.POSITIVE_INFINITY;
            ArrayList []numValues= new ArrayList[1];
            numValues[0] = new ArrayList<Double>();
            String val;
            double cval;
            for(i=0; i<totalRecord;i++)
            {
               val=datasetG[i][k];
               if(isMissing(val)==0)
               {
                    cval= Double.parseDouble(val);
                    if(cval>max)max=cval;
                    if(cval<min)min=cval;
                    numValues[0].add(cval);

               }
            }
            List allValues =GeneralFunctions.removeDuplicateValuesDouble(numValues[0]);
            double interval=GeneralFunctions.findInterval(allValues);
            double numVals = ((max-min)/interval)+1;
            int nums = (int) Math.round(numVals);

            int NofGroups=nums;
            double noI=1.0;

            NofGroups=(int)Math.sqrt((double)nums);
            int ncp=NofGroups+1;//no. of cut-points
            noI=(double)((max-min)/NofGroups);

               if(NofGroups>0)
                {
                   double diff=interval*noI;
                   double []lowDomain=new double[ncp];
                   double tmpL=min;
                   int cp=0;
                   lowDomain[cp]=tmpL;
                   tmpL=tmpL+diff;cp++;
                   for(;cp<ncp-1;cp++)
                   {
                       lowDomain[cp]=tmpL;
                       tmpL=tmpL+diff;
                   }
                   lowDomain[ncp-1]=max;


                   for(i=0; i<totalRecord;i++)
                    {
                    val=dataset[i][k];
                    if(MV[i][k]==0)
                     {
                        cval= Double.parseDouble(val);
                        int fg=-1;
                        for(j=0;j<ncp-1;j++)
                        {

                            if(j==0)
                            {
                              if(cval>=lowDomain[j] && cval<=lowDomain[j+1])
                                {
                                    fg=j;break;
                                }
                            }
                            else{
                                if(cval>lowDomain[j] && cval<=lowDomain[j+1])
                                {
                                    fg=j;break;
                                }
                            }
                         }
                        datasetG[i][k]=fg+"";
                     }
                    }
             }

         }
        }
    }



 /*
 * this method is used to generalise all numerical attributes
 * into sqrt|domainsize| categories
 */
private void generalise_old()
    {
        int i,j;

        for(int k=0; k<totalAttr;k++)
        {
         if(attributeType[k]==1)
         {
            double max=Double.NEGATIVE_INFINITY;
            double min=Double.POSITIVE_INFINITY;
            ArrayList []numValues= new ArrayList[1];
            numValues[0] = new ArrayList<Double>();
            String val;
            double cval;
            for(i=0; i<totalRecord;i++)
            {
               val=datasetG[i][k];
               if(isMissing(val)==0)
               {
                    cval= Double.parseDouble(val);
                    if(cval>max)max=cval;
                    if(cval<min)min=cval;
                    numValues[0].add(cval);

               }
            }
            List allValues =GeneralFunctions.removeDuplicateValuesDouble(numValues[0]);
            double interval=GeneralFunctions.findInterval(allValues);
            double numVals = ((max-min)/interval)+1;
            int nums = (int) Math.round(numVals);

            int NofGroups=nums;
            double noI=1.0;

            NofGroups=(int)Math.sqrt((double)nums);
            int ncp=NofGroups+1;//no. of cut-points
            noI=(double)((max-min)/NofGroups);

               if(NofGroups>0)
                {
                   double diff=interval*noI;
                   double []lowDomain=new double[ncp];
                   double tmpL=min;
                   int cp=0;
                   lowDomain[cp]=tmpL;
                   tmpL=tmpL+diff;cp++;
                   for(;cp<ncp-1;cp++)
                   {
                       lowDomain[cp]=tmpL;
                       tmpL=tmpL+diff;
                   }
                   lowDomain[ncp-1]=max;


                   for(i=0; i<totalRecord;i++)
                    {
                    val=dataset[i][k];
                    if(MV[i][k]==0)
                     {
                        cval= Double.parseDouble(val);
                        int fg=-1;
                        for(j=0;j<ncp-1;j++)
                        {

                            if(j==0)
                            {
                              if(cval>=lowDomain[j] && cval<=lowDomain[j+1])
                                {
                                    fg=j;break;
                                }
                            }
                            else{
                                if(cval>lowDomain[j] && cval<=lowDomain[j+1])
                                {
                                    fg=j;break;
                                }
                            }
                         }
                        datasetG[i][k]=fg+"";
                     }
                    }
             }

         }
        }
    }



// this method is used to print score to file
private void arrayToFile(String [][]data, int totRec, int totAttr,String outF)
{
        FileManager fileManager=new FileManager();
        fileManager.FormatDecimalPlaces(data, attributeType, MV, MR);
        File outFile=new File(outF);
        for(int i=0;i<totRec;i++)
        {
            String rec="";
            for(int j=0;j<totAttr;j++)
           {
            rec=rec+data[i][j]+", ";
           }
           if(i<totRec-1)
               rec=rec+"\n";
           if(i==0)
               fileManager.writeToFile(outFile, rec);
           else
               fileManager.appendToFile(outFile, rec);
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
private double calculateRMSE(String [][]preData,String [][]curData,int [][]mv, int []mr, int[]attrType)
{
    double rmse=0.0,c,p;
    int noR,noA;
    noA=attrType.length;
    noR=preData.length;

    double sum=0.0;
    for(int i=0;i<noR;i++)
        {
        if(mr[i]==1)
        {
            for(int j=0;j<noA;j++)
            {
                if (mv[i][j]==1)
                {
                  if (attrType[j]==1 && isMissing (preData[i][j])==0&&
                          isMissing (curData[i][j])==0)
                    {
                        p=Double.parseDouble(preData[i][j]);
                        c=Double.parseDouble(curData[i][j]);
                        sum+=Math.pow((p-c),2.0);
                     }
                  else
                  {
                      if(!preData[i][j].equals(curData[i][j]))  sum+=1.0;
                    }
                }
           }
         }
        }
    if(totalMissing>0)rmse=sum/totalMissing;
    t++;
    if(t>4)
    {
    DecimalFormat df = new DecimalFormat("####0.00000000000000000000" );
    rmse=Double.parseDouble(df.format(rmse));
    }
    return rmse;
}

 /**
  * A method to copy a 2D array
  */
 private void copy2DArray(String [][]src, String [][]dest)
 {
      for(int i=0; i<src.length;i++)
     {
         System.arraycopy(src[i], 0, dest[i], 0, src[i].length);
     }
 }
  

}
