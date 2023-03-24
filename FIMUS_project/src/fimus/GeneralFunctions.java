
package fimus;

import java.util.*;

/**
 * This class has functions that are used by several of the decision tree
 * generation methods.
 * 
 * @author helengiggins
 * @version 1.0 29/11/2010
 */
public class GeneralFunctions {

    /** Math.log(2) */
    public static final double Log2 = Math.log(2);

    /** any floating point numbers closer than one part 
     * in 10^Exponent is treated as equal
     */
    private static final int Exponent = 8;

    /** Math.pow(10.0, -Exponent) */
    public static final double Threshold = Math.pow(10.0, -Exponent);

    /** 
     * Given an array of probabilities calculates the average entropy of these
     * probabilities. Calculates the sum of <code>probs[x]</code>log<sub>2</sub>
     * <code>(1/probs[x])</code> where <code>x</code> is the number of values in
     * <code>probs</code>. Assumes all probabilities are in the range [0,1] and
     * sum to 1.0.
     * 
     * @param probs the probabilities to be used to calculate the entropy
     * @return the sum of
     * <code>probs[x]</code>log<sub>2</sub><code>(1/probs[x])</code>
     * where <code>x</code> is the number of values in <code>probs</code>
     */
    public static double calculateEntropy(double [] probs)
    {
        double entropy=0;
        /** loop over all values of probs to get the entropy */
        for(int i=0; i<probs.length; i++)
        {
            /* need to check if prob is 0 */
            if(probs[i]==0)
            {
               entropy+=0;
            }
            else
            {
                entropy+=probs[i] * (Math.log((1/probs[i])) / Log2);
            }
        }
        return entropy;
    }

    /**
     * Given an array of frequencies and a total for all frequencies, 
     * this method calculates the average entropy.
     * Calculates the sum of <code>(freq[x]/total)</code>log<sub>2</sub>
     * <code>(total/freq[x])</code> where <code>x</code> is the number of
     * values in <code>freq</code>.
     * Assumes all frequencies are in the range [0,<code>total</code>] and
     * sum to <code>total</code>.
     * 
     * @param freq an array of frequency/support counts for a set of values
     * @param total the sum of all supports
     * @return the entropy for this set of values
     */
    public static double calculateEntropy(int [] freq, int total)
    {
        double entropy=0;
        /** loop over all frequencies to get the entropy */
        int sum = 0;
        for(int i=0; i<freq.length; i++)
        {
            /* need to check if freq is 0 */
            if(freq[i]==0)
            {
               entropy+=0;
            }
            else
            {
            //note need to times denominator by 1.0 to stop an integer division
                double left =freq[i]/(total*1.0);
                double right = (total/(freq[i]*1.0));
                entropy+=(left) * (Math.log(right)/ Log2);
                sum+=freq[i];
            }
        }
 if(sum!=total)System.out.println("ERROR: Entropy calc = " + sum + " " + total);
        return entropy;
    }

    /**
     * Rounds a decimal to the specified number of decimal places.
     *
     * @param value number to be rounded
     * @param positions required decimal places
     * @return value rounded to the required number of decimal places
     */
    public static double roundToDecimals(double value, int positions)
    {
        //add a half so that when we take the floor we get the value rounded 
        //to nearest whole, rather than truncated
        double tempD = (value*Math.pow(10,positions)) + 0.5;
        
        int temp = (int)Math.floor(tempD);
        return ((double)temp)/Math.pow(10, positions);
    }

    /**
     * Given a list of {@link String} objects, sorts the list and trims 
     * any duplicate values. Note the returned list will also be ordered.
     *
     * @param unordered an unordered list of String, which may contain
     * duplicates of some values
     * @return an ordered list with only single copies of each value
     */
    public static List<String>
            removeDuplicateValuesString(ArrayList<String> unordered)
    {
  /** convert the list to an array and sort using a default java sort method */
       String [] unorderedArray = new String[unordered.size()];
       unorderedArray = unordered.toArray(unorderedArray);
       Arrays.sort(unorderedArray);
       //note: will sort in natural order for object type
       List <String> ordered = Arrays.asList(unorderedArray);
       //convert back to list
       List <String> trimmed = new ArrayList<String>();
       //add single copies of each value to this list
       /** get first value in list as current */
       String currValue = ordered.get(0);
       trimmed.add(currValue);//add first value to the list
       /** loop and remove duplicate values */
       for(int i=1; i<ordered.size(); i++)
       {
           String nextValue = ordered.get(i);
        /** check if we have a new current value, if so, store it, and add it
            * to the trimmed list
            */
           if(!currValue.equals(nextValue))
           {
               currValue = nextValue;
               trimmed.add(currValue);
           }
       }
       //for(int i=1; i<trimmed.size(); i++)
       //{
       //    System.out.println(trimmed.get(i));
       //}
       return trimmed; //return as a list
    }
  //------------
  // Constructor to generate combination
  //------------

  public long  CombinationGenerator (int n, int r)
  {
    if (n>10)
            n=4;
    else if (n>5)
            n=(int)n/2;
    else if (n>3)
            n-=(n-3);
    if (r!=2)r=2;

    long nFact = getFactorial (n);
    long rFact = getFactorial (r);
    long nminusrFact = getFactorial (n - r);
    return nFact/ (rFact*nminusrFact);
  }


 //------------------
  // Compute factorial
  //------------------

  private static long getFactorial (int n) {
    long fact = 1;
    for (int i = n; i > 1; i--) {
      fact = fact*i;
    }
    return fact;
  }


    /**
     * Given a list of {@link Double} objects, sorts the list and trims 
     * any duplicate values. Note the returned list will also be ordered.
     *
     * @param unordered an unordered list of doubles, which may contain
     * duplicates of some values
     * @return an ordered list with only single copies of each value
     */
    public static List<Double>
            removeDuplicateValuesDouble(ArrayList<Double> unordered)
    {
        /** convert the list to an array and sort using
         * a default java sort method */
       Double [] unorderedArray = new Double[unordered.size()];
       unorderedArray = unordered.toArray(unorderedArray);
       Arrays.sort(unorderedArray);
       //note: will sort in natural order for object type
       List <Double> ordered = Arrays.asList(unorderedArray);
       //convert back to list
       List <Double> trimmed = new ArrayList<Double>();
       //add single copies of each value to this list
       /** get first value in list as current */
       Double currValue = ordered.get(0);
       trimmed.add(currValue);//add first value to the list
       /** loop and remove duplicate values */
       for(int i=1; i<ordered.size(); i++)
       {
           Double nextValue = ordered.get(i);
         /** check if we have a new current value, if so, store it, and add it
            * to the trimmed list
            */
           if(!currValue.equals(nextValue))
           {
               currValue = nextValue;
               trimmed.add(currValue);
           }
       }
       //System.out.println("-o-------");
       //for(int i=0; i<ordered.size(); i++)
       //{
       //   System.out.println("ordered: " + ordered.get(i));
       //}
       //System.out.println("-t-------");
       //for(int i=0; i<trimmed.size(); i++)
       //{
       //    System.out.println(trimmed.get(i));
       //}
       return trimmed; //return as a list
    }

    /**
     * Find the minimum possible interval between two values in the list. If 
     * list size is smaller than 2, return 0 to indicate these is no interval.
     * For numbers that are really ints, find gcd(a,b) amongst all values.
     * For doubles, determine the max number of decimal places.
     *
     * @param values an ordered list of doubles containing no duplicate values
     * @return 0 when list size &lt; 2, smallest interval between
     * values otherwise
     */
    public static double findInterval(List <Double> values)
    {
        /** check we have at least two values */
        if(values.size()<2){return 0.0;} //error flag when empty or of size 1

        /** interval to 1*10^(-number of decimal places)*/
        int numDecs = getMaximumDecimalPlaces(values);
        double interval = Math.pow(10.0, numDecs*-1.0);
        //System.out.println("interval:" + interval +
        //" numDec: " + numDecs + " " + Math.pow(1, numDecs*-1.0));
        return interval;
    }

    /**
     * Determine maximum number of decimal places for any value in the list.
     *
     * @param values a list of doubles
     * @return the largest number of places after the decimal for any value
     */
    private static int getMaximumDecimalPlaces(List <Double> values)
    {
        int maxPlaces=0;
        for(int i=0; i<values.size(); i++)
        {
            double origValue = values.get(i);
            /* we can check if the current value can be converted into an int
             * when times by 10^maxPlaces. If not increase maxPlaces
             * and try again
             */
            double multiplier = Math.pow(10.0, maxPlaces);
            double currValue = origValue*multiplier;
            /** round to make sure java is not being silly */
            currValue = roundToDecimals(currValue,1);
            /** if we don't currValue as an int, add 1 to maxPlaces */
           /*
            * Original loop was:  while(Math.floor(currValue)!= currValue)
            * it becomes infinite when I used 4 decimals number for CA data set
            * 
            */

           while(Math.floor(currValue)!= currValue && maxPlaces<10)
            {
                /** NOTE: need to formulate this way due to problem with java
                 * and numbers containing 3s. Should NOT have to round
                 * currValue to 1 decimal
                 */
    //System.out.println("arggh " + Math.abs(Math.floor(currValue)-currValue));
                maxPlaces++;
                multiplier = Math.pow(10.0, maxPlaces);
                currValue = origValue*multiplier;
                currValue = roundToDecimals(currValue,1);
   //System.out.println("currValue: "+ Math.floor(currValue) + " " + currValue);
            }
        }
        return maxPlaces;
    }




}
