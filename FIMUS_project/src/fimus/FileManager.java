

package fimus;

import java.io.*;
import java.util.logging.*;
import java.util.*;


import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.channels.FileChannel;
import java.text.DecimalFormat;

/**
 * Handles all file manipulation for the application.
 *
 * @author helengiggins
 * @version 1.0 27/11/2010
 */
public class FileManager {
    final static String className = fimus.FileManager.class.getName();
    final String newline = "\n";
   

    /**
     * Class constructor.
     */
    public FileManager()
    {
        //System.out.println(className + " created");
    }

    /**
     * Reads the contents of the file to a single <code>String</code>. In the
     * event of an error the error message is returned.
     *
     * @param file the file to be read, assumes full path details
     * @return the contents of the file, unprocessed
     */
    public String readFile(File file)
    {
        String contents="";

        /** open the file and simply add each new line to the end of the
         * String, with a new line character added
         */
        try {
            FileReader fr = new FileReader(file);
            BufferedReader inFile = new BufferedReader(fr);
            String currLine = inFile.readLine();
            //keep reading until file is empty
            while(currLine!=null)
            {
                contents+=currLine + newline;
                currLine = inFile.readLine();
            }
            inFile.close();
        } catch (FileNotFoundException e) {
            String error="Error: The File " + file + " was not found, " + e;
            contents=error;
            Logger.getLogger(className).log(Level.WARNING,error);
        } catch (IOException e) {
            String error = "Error: IO Exception occured " + e;
            contents=error;
            Logger.getLogger(className).log(Level.WARNING,error);
        }
        return contents;
    }

    /**
     * Reads the contents of the file to a <code>String</code> array, where
     * each line of the file is in an array cell. In the event of an error
     * the array will be null.
     *
     * @param file file to be read, assumes full path details
     * @return each line of the file in new cell in the array
     */
    public String [] readFileAsArray(File file)
    {
        /** on reading the file each line is temporarily stored in a
         * LinkedList to determine the size of the array, then elements of the
         * LinkedList are copied to the String array.
         *
         * O(n), where n is number of lines in the file.
         */
        LinkedList <String> tempList = new LinkedList<String>();
        String [] fileArr=null;

        /** open the file and simply add each new line to the end of the list
         */
        try {
            FileReader fr = new FileReader(file);
            BufferedReader inFile = new BufferedReader(fr);
            String currLine = inFile.readLine();
            //keep reading until file is empty
            while(currLine!=null)
            {
                tempList.add(currLine);
                currLine = inFile.readLine();
            }
            inFile.close();
        } catch (FileNotFoundException e) {
            String error="Error: The File " + file + " was not found, " + e;
            Logger.getLogger(className).log(Level.WARNING,error);
        } catch (IOException e) {
            String error = "Error: IO Exception occured " + e;
            Logger.getLogger(className).log(Level.WARNING,error);
        }
        /** now copy the List elements over to the array. */
        int listSize = tempList.size();
        if(listSize>0)
        {
            fileArr = new String[listSize];
            for(int i=0; i<listSize; i++)
            {
                fileArr[i] = tempList.removeFirst();
            }
        }
        return fileArr;
    }

    /**
     * Reads the contents of the file to a <code>&lt;String&gt; List</code>,
     * each line of the file is in an element in the list. In the event of an error
     * the list will be null.
     *
     * @param file file to be read, assumes full path details
     * @return each line of the file as an element in a <code>&lt;String&gt; List</code>
     */
    public List readFileAsList(File file)
    {
        /** on reading the file each line is stored in a
         * LinkedList and returned at end of method
         *
         * O(n), where n is number of lines in the file.
         */
        List <String> tempList = new LinkedList<String>();

        /** open the file and simply add each new line to the end of the list
         */
        try {
            FileReader fr = new FileReader(file);
            BufferedReader inFile = new BufferedReader(fr);
            String currLine = inFile.readLine();
            //keep reading until file is empty
            while(currLine!=null)
            {
                tempList.add(currLine);
                currLine = inFile.readLine();
            }
            inFile.close();
        } catch (FileNotFoundException e) {
            String error="Error: The File " + file + " was not found, " + e;
            Logger.getLogger(className).log(Level.WARNING,error);
        } catch (IOException e) {
            String error = "Error: IO Exception occured " + e;
            Logger.getLogger(className).log(Level.WARNING,error);
        }
        return tempList;
    }


    /**
     * Reads the contents of the file to a <code>String</code> array, where
     * each line of the file is in an array cell. Only read the passed number
     * of lines of the file, or all of the file in the event that
     * <code>numLines</code> &gt; file size.
     * In the event of an error opening or reading the file, the array will be null.
     *
     * @param file file to be read, assumes full path details
     * @param numLines the number of lines of the file to be read
     * @return each line of the file in new cell in the array
     */
    public String [] readFileAsArray(File file, int numLines)
    {
        /** on reading the file each line is temporarily stored in a
         * LinkedList to determine the size of the array, then elements of the
         * LinkedList are copied to the String array.
         * Only read the first 'numLines' of the file.
         *
         * O(n), where n is number of lines in the file.
         */
        LinkedList <String> tempList = new LinkedList<String>();
        String [] fileArr=null;
        int currLineCount=0; //keep track of how many lines we've read so far

        /** open the file and simply add each new line to the end of the list
         * until either the end of file is reached, or we have reached 'numLines'.
         */
        try {
            FileReader fr = new FileReader(file);
            BufferedReader inFile = new BufferedReader(fr);
            String currLine = inFile.readLine();
            currLineCount++;
            //keep reading until file is empty or we've got enough lines
            while(currLine!=null && currLineCount<=numLines)
            {
                tempList.add(currLine);
                currLine = inFile.readLine();
                currLineCount++;
                //System.out.println("currLineCount: "+ currLineCount + " " + numLines);
            }
            inFile.close();
        } catch (FileNotFoundException e) {
            String error="Error: The File " + file + " was not found, " + e;
            Logger.getLogger(className).log(Level.WARNING,error);
        } catch (IOException e) {
            String error = "Error: IO Exception occured " + e;
            Logger.getLogger(className).log(Level.WARNING,error);
        }
        /** now copy the List elements over to the array. */
        int listSize = tempList.size();
        if(listSize>0)
        {
            fileArr = new String[listSize];
            for(int i=0; i<listSize; i++)
            {
                fileArr[i] = tempList.removeFirst();
            }
        }
        return fileArr;
    }

    

    /**
     * Write output <code>String</code> to the passed <code>File</code>. This
     * method will replace any current contents of the file rather than append.
     *
     * @param file output file the contents are being written to
     * @param output the new contents of the file
     * @return a message indicating success, or an error message if there is a
     * problem writing to file
     */
    public String writeToFile(File file, String output)
    {
       /** simply open a file writer and write contents to the file */
       String returnStr = "";
       try{
           FileWriter fileWriter = new FileWriter(file);
           fileWriter.write(output);
           fileWriter.flush();
           fileWriter.close();
           returnStr= "output sucessfully written to " + file.toString();
       }
       catch(IOException ex)
       {
           returnStr="Error: IO Exception occured " + ex;
           Logger.getLogger(className).log(Level.WARNING,returnStr);
       }
       return returnStr;
    }
    /**
     * Append output <code>String</code> to the passed <code>File</code>. This
     * method will retain any current contents of the file and then append the
     * output text to the end of the file.
     *
     * @param file output file the contents are being written to
     * @param output the new contents to be added to the file
     * @return a message indicating success, or an error message if there is a
     * problem writing to file
     */
    public String appendToFile(File file, String output)
    {
       /** simply open a file writer and write contents to the file */
       String returnStr = "";
       try{
           FileWriter fileWriter = new FileWriter(file,true);
           fileWriter.write(output);
           fileWriter.flush();
           fileWriter.close();
           returnStr= "output sucessfully written to " + file.toString();
       }
       catch(IOException ex)
       {
           returnStr="Error: IO Exception occured " + ex;
           Logger.getLogger(className).log(Level.WARNING,returnStr);
       }
       return returnStr;
    }

    /**
     * This method will read a file and remove any records (lines) containing
     * missing values. Missing values are denoted using a '?'. The method also
     * writes the records to the output file.
     *
     * @param file name of input file
     * @param output name of output file
     */
    public void removeMissingValuesFromFile(File file, String output)
    {
        /** on reading the file each line is stored in a
         * LinkedList and returned at end of method. Records are shuffled before
         * being returned, so they are in a random order.
         *
         * Note, any line read containing a '?' is ignored. The question mark is
         * used to indicate a missing value.
         *
         * O(n), where n is number of lines in the file.
         */
        List <String> tempList = new LinkedList<String>();

        /** open the file and simply add each new line to the end of the list
         */
        try {
            FileReader fr = new FileReader(file);
            BufferedReader inFile = new BufferedReader(fr);
            String currLine = inFile.readLine();
            int totalLines = 0;
            int removedLines = 0;
            //keep reading until file is empty
            while(currLine!=null)
            {
                /** test for missing value, only add line if none is present */
                if(currLine.indexOf("?")<0)
                {
                    tempList.add(currLine);
                    
                }
                else
                {
                    removedLines++;
                    
                }
                currLine = inFile.readLine();
                totalLines++;
            }
            inFile.close();
//            System.out.println("Lines in original file: " + totalLines + "\n"+
//                    "Lines removed: " + removedLines + "\n");//commentd by gea

        } catch (FileNotFoundException e) {
            String error="Error: The File " + file + " was not found, " + e;
            Logger.getLogger(className).log(Level.WARNING,error);
        } catch (IOException e) {
            String error = "Error: IO Exception occured " + e;
            Logger.getLogger(className).log(Level.WARNING,error);
        }
        /** extract Strings and then write to file */
        StringBuilder outStr = new StringBuilder();
        for(int i=0; i<tempList.size(); i++)
        {
          outStr.append(tempList.get(i));
          outStr.append("\n");
        }
        //write the records to file
        this.writeToFile(new File(output), outStr.toString());
    }

    /**
     * Given a data file and number of folds for cross validation, generate that
     * number of training-testing data set pairs.
     *
     * @param file the input file name
     * @param out the output file name (number will be appended to this to distinguish files)
     * @param num the number of folds for cross-validation
     */
    public void generateCrossFoldValidationFiles(File file, File out, int num)
    {
        /** read in the file and store to a list */
        List <String> originalFile = readFileAsList(file);
        /** shuffle the list so that records are in no particular order */
        Collections.shuffle(originalFile);
        /** determine the number of records in each fold */
        int totalRecords = originalFile.size();
        int foldNum = totalRecords/num;
        int lastFoldNum = foldNum + (totalRecords - (foldNum*num));
//        System.out.println("TotalRecords:" + totalRecords + " Fold 1-(num-1):" + foldNum + " Fold num: " + lastFoldNum);
        int currStartPos = 0; //will increment after each of the set of records is extracted

        /** now loop over the number of folds to be generated
         * creating two lists, one for training data and one for testing
         *
         */
        for(int i=0; i<num; i++)
        {
            ArrayList<String> training = new ArrayList<String>();
            ArrayList<String> testing = new ArrayList<String>();
            for(int j=0;j<totalRecords; j++)
            {
                /** if not at last fold, which can be a different size */
                if(i!=(num-1))
                {
                    /** if we're in the right range for test data */
                    if((j>=currStartPos) && (j<(currStartPos+foldNum)))
                    {
                        //System.out.println(i + " " + j + " X");
                        testing.add(originalFile.get(j));
                    }
                    else /*training */
                    {
                        //System.out.println(i + " " + j + " O");
                        training.add(originalFile.get(j));
                    }
                }
                else
                {
                    /** if we're in the right range for test data */
                    if((j>=currStartPos) && (j<(currStartPos+lastFoldNum)))
                    {
                        //System.out.println(i + " " + j + " X");
                        testing.add(originalFile.get(j));
                    }
                    else /*training */
                    {
                        //System.out.println(i + " " + j + " O");
                        training.add(originalFile.get(j));
                    }
                }
            }
            /**
             * Write the training and testing lists to file
             */
            /** create a new file with '-x' in front of output file name, where x is the fold number */
            File outFile = out;
            String outFilePath = outFile.getPath();
            int indexOfDot = outFilePath.lastIndexOf("."); //to get position of file extension
            String pathBeforeExtension = outFilePath.substring(0,indexOfDot);
            String trainOutName = pathBeforeExtension + "-train-" + (i+1) + ".txt";
            outFile = new File(trainOutName);
            /** extract Strings and then write to file */
            StringBuilder outStr = new StringBuilder();
            for(int x=0; x<training.size(); x++)
            {
              outStr.append(training.get(x));
              outStr.append("\n");
            }
            //write the records to file
            this.writeToFile(outFile, outStr.toString());
            /** create a new file with '-x' in front of output file name, where x is the fold number */
            outFile = out;
            outFilePath = outFile.getPath();
            indexOfDot = outFilePath.lastIndexOf("."); //to get position of file extension
            pathBeforeExtension = outFilePath.substring(0,indexOfDot);
            String testOutName = pathBeforeExtension + "-test-" + (i+1) + ".txt";
            outFile = new File(testOutName);
            /** extract Strings and then write to file */
            outStr = new StringBuilder();
            for(int x=0; x<testing.size(); x++)
            {
              outStr.append(testing.get(x));
              outStr.append("\n");
            }
            //write the records to file
            this.writeToFile(outFile, outStr.toString());
            currStartPos+=foldNum; //update the startPosition for the testing set
        }
    }

    /**
     * Generate a name file from a data file. Also need to input a two line file
     * containing the attribute types one first line, and attribute names on second.
     * Line 1 shows attribute type for each attribute,
     * 0=categorical, 1=numerical, 2=class e.g.
     * <p>l 0 1 2 0 1</p>
     * <p>Line 2 shows attribute names, note that names can be separated by spaces,
     * tabs, or commas. However, names need to be one word, and can contain any
     * character except for comma.</p>
     * @param attrInfo the file containing attribute info on two lines
     * @param dataFile the data file we are generating the name file for
     * @param outFile the filename for the new name file
     * @return an message, either indicating success, or describing error
     */
    public String extractNameFileFromDataFile(File attrInfo, File dataFile, File outFile)
    {
        /** read the two files and store as arrays */
        String [] attrFile = readFileAsArray(attrInfo, 2);
        String [] dataStrings = readFileAsArray(dataFile);
        /** store return message */
        String retStr="Resulting name file successfully written to " + outFile.getPath();
        /** read the first line of the attribute info file to determine the number
         * and type of attributes, as well as class attribute
         */
        int [] attrTypes=null; //will store attribute type for each attribute
        String [] attrNames = null; //store the names of each attribute
        int classIndex = -1; //the attribute index of the class attribute
        int numAttrs = -1; //the number of attributes
        double [] numHighDomain=null; //store highest value for numerical attrs
        double [] numLowDomain=null;  //store lowest value for numerical attrs
        double [] intervals = null; //store the minimum interval between values
        ArrayList [] catValues = null; //for each categorical attr store a list of string values
        ArrayList [] numValues = null; //store list of numerical values
        try{
            /** tokenize first line to find out how many attributes we have, and then
             * store their type, tokenize second line to get attribute names
             */
            StringTokenizer tokens = new StringTokenizer(attrFile[0], " ,\n\t"); //remove spaces, commas, tabs and newlines
            StringTokenizer tokensNames =new StringTokenizer(attrFile[1], " ,\n\t"); //remove spaces, commas, tabs and newlines
            numAttrs = tokens.countTokens();
            attrTypes = new int [numAttrs];
            attrNames = new String[numAttrs];
            /** initalize variables for storing value info */
            numHighDomain = new double[numAttrs];
            numLowDomain = new double[numAttrs];
            intervals = new double[numAttrs];
            catValues = new ArrayList[numAttrs];
            numValues = new ArrayList[numAttrs];
            /** read attribute info, and initialize appropriate variables*/
            for(int currAttr=0; currAttr<numAttrs; currAttr++)
            {
                int currAttrType = Integer.parseInt(tokens.nextToken());//get attr type
                attrTypes[currAttr]=currAttrType;
                attrNames[currAttr]=tokensNames.nextToken(); //get attr name
                /** check for attribute type */
                if(currAttrType==2)//class attribute
                {
                    classIndex = currAttr;
                    catValues[currAttr] = new ArrayList<String>();
                }
                else if(currAttrType==0)//categorical attribute
                {
                    catValues[currAttr] = new ArrayList<String>();
                }
                else //numerical attribute
                {
                    numValues[currAttr] = new ArrayList<Double>();
                    numHighDomain[currAttr] = Double.NEGATIVE_INFINITY;//defaults
                    numLowDomain[currAttr] = Double.POSITIVE_INFINITY;//defaults
                }
            }
        }
        catch(Exception e)
        {
            retStr="There was a problem extracting data from the attribute info file.\n"+
                    "Please check the file and try again, no name file created.\n"+ e.toString();
            return retStr;
        }
        /** 
         * Now we know how many attributes and types, we need to find the ranges and values.
         * For categorical attributes, we are building a list of values. For numerical
         * values, we are trying to find highest and lowest value, and storing all values
         * to an Integer list. Later we will sort these list to determine the minimum interval
         * between values, and the values for categorical.
         */
        for(int currRec=0; currRec<dataStrings.length; currRec++)
        {
            //System.out.println("rec: " + currRec + "dataStrings " + dataStrings.length + "\n" + dataStrings[currRec]);
            /**for each record, we need to tokenize the values for each attribute */

            /** NOTE: Comment in this 'if' statement to read only every second line of the file **/
            //if(currRec%2==0){
             StringTokenizer tokens = new StringTokenizer(dataStrings[currRec], " ,\t\n\r\f");
             //this.appendToFile(new File("../centNervNew.txt"),dataStrings[currRec]+"\n" ); //used to output every second line of file
            /** for each attribute value, we determine which type */
            for(int currAttr=0; currAttr<numAttrs; currAttr++)
            {       
                String currValue = tokens.nextToken();
                //System.out.println("currAttr: " + currAttr + " currValue: " + currValue);
                /** if categorical attribute, or class attribute, just add the value
                 * to the list for that attribute
                 */
                if(attrTypes[currAttr]==0 || attrTypes[currAttr]==2)
                {
                    catValues[currAttr].add(currValue);
                }
                /** if numerical, convert the String value to a double, and store to list.
                 * Also check if we have a new high or low value.
                 */
                else if(attrTypes[currAttr]==1)
                {
                    double numValue = 0.0;
                    if(isMissing(currValue)==0)
                            numValue =Double.parseDouble(currValue);
                    double currMax = numHighDomain[currAttr];
                    double currLow = numLowDomain[currAttr];
                    //System.out.println(numLowDomain[currAttr] + " " + numHighDomain[currAttr]);
                    if(numValue>currMax)
                    {
                        numHighDomain[currAttr]=numValue;
                    }
                    if(numValue<currLow)
                    {
                        numLowDomain[currAttr]=numValue;
                    }
                    numValues[currAttr].add(numValue);
                }

            }
            //}// uncomment for if statement to appear
        }

       /** remove duplicate values from each list, and find intervals for numerical
        * values
        */
        List [] allValues = new List[numAttrs];
        for(int currAttr=0; currAttr<numAttrs; currAttr++)
        {
            
            /** if categorical attribute, or class attribute
             */
            if(attrTypes[currAttr]==0 || attrTypes[currAttr]==2)
            {
                allValues[currAttr] = GeneralFunctions.removeDuplicateValuesString(catValues[currAttr]);
            }
            /** if numerical, convert the String value to a double, and store to list.
             * Also check if we have a new high or low value.
             */
            else if(attrTypes[currAttr]==1)
            {
                allValues[currAttr] = GeneralFunctions.removeDuplicateValuesDouble(numValues[currAttr]);
                /** find intervals */
                intervals[currAttr] = GeneralFunctions.findInterval(allValues[currAttr]);
            }
        }
        /** Now have all info we need to name file, just need to build the output String */
        StringBuilder outStr = new StringBuilder();
        /*<p>Note on format of nameFile</p>
         * <ul>
         *   <li><strong>First line:</strong> class attribute index, number of class values</li>
         *   <li><strong>Second line:</strong> number of records, number of attributes</li>
         *   <li><strong>Categorical attribute:</strong> <code>c</code>, attribute name, number of categories, values</li>
         *   <li><strong>Numerical attribute:</strong> <code>n</code>, attribute name, low domain,
         *       high domain, interval, number of values</li>
         * </ul>
         */
        /** first line */
        int numClasses = allValues[classIndex].size();
        String firstLine = classIndex + ", " + numClasses + ",\n";
        outStr.append(firstLine);
        /** second line */
        String secondLine = dataStrings.length + ", " + numAttrs + ",\n";
        outStr.append(secondLine);
        /** now for each attribute */
        for(int currAttr=0; currAttr<numAttrs; currAttr++)
        {
            StringBuilder currLine= new StringBuilder(); //better for longer strings
            /** categorical or class attribute */
            if(attrTypes[currAttr]==0 || attrTypes[currAttr]==2)
            {
                currLine.append("c, ");
                currLine.append(attrNames[currAttr]); //attr name
                currLine.append(", ");
                List currCatAttr = allValues[currAttr];
                int numCats = currCatAttr.size();
                currLine.append(numCats); //number of categories
                currLine.append(", ");
                /** now print values */
                for(int i=0; i<numCats; i++)
                {
                    currLine.append(currCatAttr.get(i));
                    currLine.append(", ");
                }

            }
            /* numerical */
            else
            {
                currLine.append("n, ");
                currLine.append(attrNames[currAttr]); //attr name
                currLine.append(", ");
                double numVals = ((numHighDomain[currAttr]-numLowDomain[currAttr])/intervals[currAttr])+1;
                int nums = (int) Math.round(numVals);
                /*low domain, high domain, interval, number of values */
                String otherDetails = numLowDomain[currAttr]+ ", "+numHighDomain[currAttr]+ ", "
                        + intervals[currAttr]+ ", " + nums + ",";
                currLine.append(otherDetails);
            }
            currLine.append("\n");//end line
            outStr.append(currLine.toString());
        }
        //write the output string to file
        this.writeToFile(outFile, outStr.toString());

        return retStr;
    }

    /*
     * the following three methods has been written by gea
     */
    /**
     * This method will read a file and remove any records (lines) containing
     * missing values. Missing values are denoted using a '?'. The method also
     * writes the records to the output file.
     *
     * @param file name of input file
     * @param output_NM name of output file having NO missing values
     * @param output_MM name of output file having missing values
     */
public void divideDataset(File file, String output_NM,String output_MM)
    {
        /** on reading the file each line is stored in a
         * LinkedList and returned at end of method. Records are shuffled before
         * being returned, so they are in a random order.
         *
         * Note, any line read containing a '?' is ignored. The question mark is
         * used to indicate a missing value.
         *
         * O(n), where n is number of lines in the file.
         */
        List <String> tempList = new LinkedList<String>();
        List <String> tempList_MM = new LinkedList<String>();
        /** open the file and simply add each new line to the end of the list
         */
        try {
            FileReader fr = new FileReader(file);
            BufferedReader inFile = new BufferedReader(fr);
            String currLine = inFile.readLine();
            int totalLines = 0;
            int removedLines = 0;
            //keep reading until file is empty
            while(currLine!=null)
            {
                /** test for missing value, only add line if none is present */
                if(currLine.indexOf("?")<0)
                {
                    tempList.add(currLine);

                }
                else
                {
                    tempList_MM.add(currLine);
                    removedLines++;

                }
                currLine = inFile.readLine();
                totalLines++;
            }
            inFile.close();
//            System.out.println("Lines in original file: " + totalLines + "\n"+
//                    "Lines removed: " + removedLines + "\n");

        } catch (FileNotFoundException e) {
            String error="Error: The File " + file + " was not found, " + e;
            Logger.getLogger(className).log(Level.WARNING,error);
        } catch (IOException e) {
            String error = "Error: IO Exception occured " + e;
            Logger.getLogger(className).log(Level.WARNING,error);
        }
        /** extract Strings and then write to file */
        StringBuilder outStr = new StringBuilder();
        for(int i=0; i<tempList.size(); i++)
        {
          outStr.append(tempList.get(i));
          outStr.append("\n");
        }
        //write the records to file
        this.writeToFile(new File(output_NM), outStr.toString());
         /** extract Strings and then write to file */
        outStr = new StringBuilder();
        for(int i=0; i<tempList_MM.size(); i++)
        {
          outStr.append(tempList_MM.get(i));
          outStr.append("\n");
        }
        //write the records to file
        this.writeToFile(new File(output_MM), outStr.toString());
    }
/**
     * This method will read a file and write it to the output file.
     *
     * @param srFile name of source file
     * @param dtFile name of output file
     */
public void copyFile(String srFile, String dtFile){
    try{
      File f1 = new File(srFile);
      File f2 = new File(dtFile);
      InputStream in = new FileInputStream(f1);
     //For Overwrite the file.
      OutputStream out = new FileOutputStream(f2);

      byte[] buf = new byte[1024];
      int len;
      while ((len = in.read(buf)) > 0){
        out.write(buf, 0, len);
      }
      in.close();
      out.close();

    }
    catch(FileNotFoundException ex){
      System.out.println(ex.getMessage() + " in the specified directory.");
      System.exit(0);
    }
    catch(IOException e){
      System.out.println(e.getMessage());
    }
  }

/**
     * Reads the contents of the file to a <code>String</code> array, where
     * each line of the file is in an array cell. In the event of an error
     * the array will be null.
     *
     * @param file file to be read, assumes full path details
     * @return each attribute value of the file in new cell in the array
     */
    public void FormatDecimalPlaces(String [][]data, int[]attrNtype,
            int[][]MV,int[]MR)
    {
        int tr=data.length;int tm=0;
        int attr=attrNtype.length;
        
        int []mdp=attrMaxDecimalPlaces(attrNtype,data);
        String []fStr=new String[attr];
        for (int c=0;c<attr;c++)
        {
            fStr[c]="";
            if (attrNtype[c]==1)
            {
                String fs="";
                for(int i=0;i<mdp[c];i++)
                {
                    fs=fs+"0";
                }
                if (mdp[c]>0)
                {
                    fStr[c] = "####0." + fs;}
                else
                {
                    fStr[c] = "####0";}
            }
        }
        for(int i=0;i<tr;i++)
        {
            for(int j=0;j<attr;j++)
            {
                if(MV[i][j]==1)tm++;
                
            }
        }
        
         for(int i=0;i<tr;i++)
        {
            for(int j=0;j<attr;j++)
            {
                if(MV[i][j]==1&& attrNtype[j]==1&& isMissing(data[i][j])==0)
                {
                     DecimalFormat df = new DecimalFormat(fStr[j]);
                     data[i][j]=df.format(Double.parseDouble(data[i][j]))+"";
                }
            }
        }


    }

/**
     * This method will be used to change a file name by a supplied padding .

     *Author Geaur
     * @param filename file name to be changed
     * @param padding to add at the right of the filename
     * @return a changed file name, or an error message if there is a
     * problem writing to file
     */
    public String changedFileName(String filename, String padding)
    {
        String returnStr = "";
       /** simply rename a filename*/
        File outFile = new File(filename);
        String outFilePath = outFile.getPath();
        int indexOfDot = outFilePath.lastIndexOf("."); //to get position of file extension
        String pathBeforeExtension = outFilePath.substring(0,indexOfDot);
        String pathAfterExtension = outFilePath.substring(indexOfDot+1,outFilePath.length());
        returnStr = pathBeforeExtension + padding+"."+pathAfterExtension;
        return returnStr;
    }

    public void removeTMPFiles(String extension)
    {
        boolean bool=false;
        java.io.File currentDir = new java.io.File("");
        File dir = new File(currentDir.getAbsolutePath());
        File[] files = dir.listFiles();
        for(int i=0; i<files.length;i++)
        {

            if(files[i].isFile())
            {
                String outFilePath = files[i].getPath();
                int indexOfDot = outFilePath.lastIndexOf("."); //to get position of file extension
                String pathAfterExtension = outFilePath.substring(indexOfDot+1,outFilePath.length());
                if(extension.equals(pathAfterExtension))
                {
                    if(files[i].exists())
                    {
                       bool= files[i].delete();
                    }
                }
            }
        }
    }
    /**
     * Reads the contents of the file to a <code>String</code> array, where
     * each line of the file is in an array cell. In the event of an error
     * the array will be null.
     *
     * @param file file to be read, assumes full path details
     * @return each attribute value of the file in new cell in the array
     */
    public String [][] readFileAs2DArray(File file)
    {
        /** on reading the file each line is temporarily stored in a
         * LinkedList to determine the size of the array, then elements of the
         * LinkedList are copied to the String array.
         *
         * O(n), where n is number of lines in the file.
         */
        LinkedList <String> tempList = new LinkedList<String>();
        String [][] fileArr=null;
        int totalAttr=0;
        StringTokenizer tokenizer;
        String currLine;
        /** open the file and simply add each new line to the end of the list
         */
        try {
            FileReader fr = new FileReader(file);
            BufferedReader inFile = new BufferedReader(fr);
            currLine = inFile.readLine();
            tokenizer= new StringTokenizer(currLine, " ,\t\n\r\f");
            totalAttr=tokenizer.countTokens();
            //keep reading until file is empty
            while(currLine!=null)
            {
                tempList.add(currLine);
                currLine = inFile.readLine();
            }
            inFile.close();
        } catch (FileNotFoundException e) {
            String error="Error: The File " + file + " was not found, " + e;
            Logger.getLogger(className).log(Level.WARNING,error);
        } catch (IOException e) {
            String error = "Error: IO Exception occured " + e;
            Logger.getLogger(className).log(Level.WARNING,error);
        }
        /** now copy the List elements over to the array. */
        int listSize = tempList.size();
        if(listSize>0)
        {
            fileArr = new String[listSize][totalAttr];
            for(int i=0; i<listSize; i++)
            {
                currLine=tempList.removeFirst();
                tokenizer= new StringTokenizer(currLine, " ,\t\n\r\f");
                for(int j=0;j<totalAttr;j++)
                    fileArr[i][j] = tokenizer.nextToken();
            }
        }
        return fileArr;
    }
    
    /**
     * Find the maximum decimal places of each attribute
     *
     * @param attrF int 1D array of the attribute info
     * @param dataF string 2D array of the data file
     * @return fStr contains max decimal places of each attribute
     */
    public int []attrMaxDecimalPlaces(int []attrF,String [][]dataF)
    {
        int totRow=dataF.length;
        int totAttr=attrF.length;
        int []fStr=new int[totAttr];
        for(int c=0; c<totAttr;c++)
        {
            int md=0;
            if(attrF[c]==1)
            {
                for (int r=0; r<totRow;r++)
                {
                    String s = dataF[r][c];
                    int indexOfDot = s.lastIndexOf("."); //to get position of file extension
                    if( indexOfDot>-1)
                    {
                    String tmps = s.substring(indexOfDot+1,s.length());
                    int l=tmps.length();
                    if (l>md)
                        md=l;
                    }
                }
            }
            fStr[c]=md;
        }
        return fStr;
    }

    /**
     * this function will remove all temporarily created files
     * 
     * @param fileName the name of the file to be removed
     */
    public void removeFile(String fileName)
    {
       File tFile;boolean bool=false;
       
        tFile=new File(fileName);
        if(tFile.exists())
        {
           bool= tFile.delete();
        }
    }
     /**
     * this function will remove all temporarily created files
     *
     * @param fileNames the names of the files to be removed
     * @param noOfFiles the number of files to be removed
     */
    public void removeListOfFiles(String []fileNames, int noOfFiles)
    {
       File tFile;boolean bool=false;
       for(int i=0; i<noOfFiles;i++)
       {
        tFile=new File(fileNames[i]);
        if(tFile.exists())
        {
           bool= tFile.delete();
        }
       }
    }

    /*
     * Convert a string into an Array
     */
public String[] StrArray(String str)
{
    StringTokenizer testR= new StringTokenizer(str, " ,\t\n\r\f");
    int totAttr=testR.countTokens();
    String []rec=new String[totAttr];
    for(int k=0;k<totAttr;k++)
    {
        rec[k]=testR.nextToken();
    }
    return rec;
}
/*
     * Convert a string into an Array
     */
public String ArrayStr(String []Ary)
{
    int totAttr=Ary.length;
    String rec="";
    for(int k=0;k<totAttr;k++)
    {
        rec=rec+Ary[k]+", ";
    }
    return rec;
}
 /*
  * this function will indicate whether or not a value is missing.
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
 
//
// The destination directory to copy to. This directory
// doesn't exists and will be created during the copy
// directory process.
 public void copyDirectory(String source, String destination)
 {
    //
        // An existing directory to copy.
        //
        File srcDir = new File(source);

        //
        // The destination directory to copy to. This directory
        // doesn't exists and will be created during the copy
        // directory process.
        //
        File destDir = new File(destination);

        try {
            //
            // Copy source directory into destination directory
            // including its child directories and files. When
            // the destination directory is not exists it will
            // be created. This copy process also preserve the
            // date information of the file.
            //
            CopyDirectoryFiles(srcDir, destDir);
        } catch (IOException e) {
            
        }
 }

 public void CopyDirectoryFiles(File sourceLocation , File targetLocation)
    throws IOException {

        if (sourceLocation.isDirectory()) {
            if (!targetLocation.exists()) {
                targetLocation.mkdir();
            }

            String[] children = sourceLocation.list();
            for (int i=0; i<children.length; i++) {
                CopyDirectoryFiles(new File(sourceLocation, children[i]),
                        new File(targetLocation, children[i]));
            }
        } else {

            InputStream in = new FileInputStream(sourceLocation);
            OutputStream out = new FileOutputStream(targetLocation);

            // Copy the bits from instream to outstream
            byte[] buf = new byte[1024];
            int len;
            while ((len = in.read(buf)) > 0) {
                out.write(buf, 0, len);
            }
            in.close();
            out.close();
        }
    }

 public void copyFile(File sourceFile, File destinationFile) {
               
                try {
                        FileInputStream fileInputStream = new FileInputStream(sourceFile);
                        FileOutputStream fileOutputStream = new FileOutputStream(
                                        destinationFile);

                        int bufferSize;
                        byte[] bufffer = new byte[512];
                        while ((bufferSize = fileInputStream.read(bufffer)) > 0) {
                                fileOutputStream.write(bufffer, 0, bufferSize);
                        }
                        fileInputStream.close();
                        fileOutputStream.close();
                } catch (Exception e) {
                    System.out.println("Error in copying file..");
                       //do nothing
                }
        }

 public void copyFile1(File sourceFile, File destFile) throws IOException {
  if (!destFile.exists()) {
    destFile.createNewFile();
  }

  FileInputStream fIn = null;
  FileOutputStream fOut = null;
  FileChannel source = null;
  FileChannel destination = null;
  try {
    fIn = new FileInputStream(sourceFile);
    source = fIn.getChannel();
    fOut = new FileOutputStream(destFile);
    destination = fOut.getChannel();
    long transfered = 0;
    long bytes = source.size();
    while (transfered < bytes) {
      transfered += destination.transferFrom(source, 0, source.size());
      destination.position(transfered);
    }
  } finally {
    if (source != null) {
      source.close();
    } else if (fIn != null) {
      fIn.close();
    }
    if (destination != null) {
      destination.close();
    } else if (fOut != null) {
      fOut.close();
    }
  }
}

}
