//  Distance.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.comparators.ObjectiveComparator;
import jmetal.util.wrapper.XReal;

/**
 * This class implements some utilities for calculating distances
 */
public class Distance {

  /**
   * Constructor.
   */
  public Distance() {
    //do nothing.
  } // Distance


  /**
   * Returns a matrix with distances between solutions in a
   * <code>SolutionSet</code>.
   * @param solutionSet The <code>SolutionSet</code>.
   * @return a matrix with distances.
   */
  public double [][] distanceMatrix(SolutionSet solutionSet) {
    Solution solutionI, solutionJ;

    //The matrix of distances
    double [][] distance = new double [solutionSet.size()][solutionSet.size()];
    //-> Calculate the distances
    for (int i = 0; i < solutionSet.size(); i++){
      distance[i][i] = 0.0;
      solutionI = solutionSet.get(i);
      for (int j = i + 1; j < solutionSet.size(); j++){
        solutionJ = solutionSet.get(j);
        distance[i][j] = this.distanceBetweenObjectives(solutionI,solutionJ);
        distance[j][i] = distance[i][j];
      } // for
    } // for        

    //->Return the matrix of distances
    return distance;
  } // distanceMatrix

  /** Returns the minimum distance from a <code>Solution</code> to a
   * <code>SolutionSet according to the objective values</code>.
   * @param solution The <code>Solution</code>.
   * @param solutionSet The <code>SolutionSet</code>.
   * @return The minimum distance between solution and the set.
   * @throws JMException
   */
  public double distanceToSolutionSetInObjectiveSpace(Solution    solution,
                                                      SolutionSet solutionSet) throws JMException{
    //At start point the distance is the max
    double distance = Double.MAX_VALUE;

    // found the min distance respect to population
    for (int i = 0; i < solutionSet.size();i++){
      double aux = this.distanceBetweenObjectives(solution,solutionSet.get(i));
      if (aux < distance)
        distance = aux;
    } // for

    //->Return the best distance
    return distance;
  } // distanceToSolutionSetinObjectiveSpace

  /** Returns the minimum distance from a <code>Solution</code> to a 
   * <code>SolutionSet according to the encodings.variable values</code>.
   * @param solution The <code>Solution</code>.
   * @param solutionSet The <code>SolutionSet</code>.
   * @return The minimum distance between solution and the set.
   * @throws JMException
   */
  public double distanceToSolutionSetInSolutionSpace(Solution    solution,
                                                     SolutionSet solutionSet) throws JMException{
    //At start point the distance is the max
    double distance = Double.MAX_VALUE;

    // found the min distance respect to population
    for (int i = 0; i < solutionSet.size();i++){
      double aux = this.distanceBetweenSolutions(solution,solutionSet.get(i));
      if (aux < distance)
        distance = aux;
    } // for

    //->Return the best distance
    return distance;
  } // distanceToSolutionSetInSolutionSpace

  /** Returns the distance between two solutions in the search space.
   *  @param solutionI The first <code>Solution</code>.
   *  @param solutionJ The second <code>Solution</code>.
   *  @return the distance between solutions.
   * @throws JMException
   */
  public double distanceBetweenSolutions(Solution solutionI, Solution solutionJ)
          throws JMException{
    /*
    double distance = 0.0;
    if ((solutionI.getDecisionVariables() != null) &&
            (solutionJ.getDecisionVariables() != null)) {
      Variable[] decisionVariableI = solutionI.getDecisionVariables();
      Variable[] decisionVariableJ = solutionJ.getDecisionVariables();

      double diff;    //Auxiliar var
      //-> Calculate the Euclidean distance
      for (int i = 0; i < decisionVariableI.length; i++){
        diff = decisionVariableI[i].getValue() -
                decisionVariableJ[i].getValue();
        distance += Math.pow(diff,2.0);
      } // for    
    }
    //-> Return the euclidean distance
    return Math.sqrt(distance);
    */
    double distance = 0.0;
    XReal solI = new XReal(solutionI) ;
    XReal solJ = new XReal(solutionJ) ;

      double diff;    //Auxiliar var
      //-> Calculate the Euclidean distance
      for (int i = 0; i < solI.getNumberOfDecisionVariables(); i++){
        diff = solI.getValue(i) - solJ.getValue(i);
        distance += Math.pow(diff,2.0);
      } // for
    //-> Return the euclidean distance
    return Math.sqrt(distance);
  } // distanceBetweenSolutions

  /** Returns the distance between two solutions in objective space.
   *  @param solutionI The first <code>Solution</code>.
   *  @param solutionJ The second <code>Solution</code>.
   *  @return the distance between solutions in objective space.
   */
  public double distanceBetweenObjectives(Solution solutionI, Solution solutionJ){
    double diff;    //Auxiliar var
    double distance = 0.0;
    //-> Calculate the euclidean distance
    for (int nObj = 0; nObj < solutionI.getNumberOfObjectives();nObj++){
      diff = solutionI.getObjective(nObj) - solutionJ.getObjective(nObj);
      distance += Math.pow(diff,2.0);
    } // for   

    //Return the euclidean distance
    return Math.sqrt(distance);
  } // distanceBetweenObjectives.

  /**
   * Return the index of the nearest solution in the solution set to a given solution
   * @param solution
   * @param solutionSet
   * @return  The index of the nearest solution; -1 if the solutionSet is empty
   */
  public int indexToNearestSolutionInSolutionSpace(Solution solution, SolutionSet solutionSet) {
    int index = -1 ;
    double minimumDistance = Double.MAX_VALUE ;
    try {
      for (int i = 0 ; i < solutionSet.size(); i++) {
        double distance = 0;
          distance = distanceBetweenSolutions(solution, solutionSet.get(i));
          if (distance < minimumDistance) {
            minimumDistance = distance ;
            index = i ;
        }
      }
    }
    catch (Exception e) {
      e.printStackTrace();
    }
    return index ;
  }

  /** Assigns crowding distances to all solutions in a <code>SolutionSet</code>.
   * @param solutionSet The <code>SolutionSet</code>.
   * @param nObjs Number of objectives.
   */
  public void crowdingDistanceAssignment(SolutionSet solutionSet, int nObjs) {
    int size = solutionSet.size();

    if (size == 0)
      return;

    if (size == 1) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if

    if (size == 2) {
      solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      solutionSet.get(1).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    } // if       

    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionSet.get(i));
    }

    for (int i = 0; i < size; i++)
      front.get(i).setCrowdingDistance(0.0);

    for (int i = 0; i<nObjs; i++) {
      // Sort the population by Obj n            
      front.sort(new ObjectiveComparator(i));
      double objetiveMinn = front.get(0).getObjective(i);
      double objetiveMaxn = front.get(front.size()-1).getObjective(i);

      //Set de crowding distance            
      front.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      front.get(size-1).setCrowdingDistance(Double.POSITIVE_INFINITY);

      for (int j = 1; j < size-1; j++) {
        double distance = front.get(j+1).getObjective(i) - front.get(j-1).getObjective(i);
        distance /= (objetiveMaxn - objetiveMinn);
        distance += front.get(j).getCrowdingDistance();
        front.get(j).setCrowdingDistance(distance);
      } // for
    } // for        
  } // crowdingDistanceAssing            

  public void FuzzyDistanceAssignment(SolutionSet solutionList, int nObjs) {
    int size = solutionList.size();
    int numberOfObjectives = nObjs;
    if (size <= nObjs) {
      for(int x=0;x<size;x++)
      solutionList.get(x).setCrowdingDistance(Double.POSITIVE_INFINITY);
      return;
    }




    //Use a new SolutionSet to evite alter original solutionSet
    SolutionSet front = new SolutionSet(size);
    for (int i = 0; i < size; i++){
      front.add(solutionList.get(i));
    }

    for (int i = 0; i < size; i++) {
      front.get(i).setCrowdingDistance(0.0);
    }

    
    double objetiveMaxn[]=new double[numberOfObjectives];
    double objetiveMinn[]=new double[numberOfObjectives];
    double distan;

    for (int i = 0; i < numberOfObjectives; i++) {
      // Sort the population by Obj n
      front.sort(new ObjectiveComparator(i));
      objetiveMinn[i] = front.get(0).getObjective(i);
      objetiveMaxn[i] = front.get(front.size() - 1).getObjective(i);

      //Set de crowding distance Los extremos si infinitos
      front.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
      front.get(size - 1).setCrowdingDistance(Double.POSITIVE_INFINITY);
      
    }
    double [][] distance = distanceMatrix_normalize(front,objetiveMaxn,objetiveMinn);
    double dminn,dmaxx;
    dminn=Double.MAX_VALUE;
    dmaxx=0.0;
    for(int i=0;i<distance.length;i++)
        for(int j=0;j<distance.length;j++)
        {
            if(i!=j){
            if(distance[i][j]<dminn)//dmin[i]
            {
                dminn=distance[i][j];
            }
            if(distance[i][j]>dmaxx)//dmax[i]
            {
                dmaxx=distance[i][j];
            }
            }
            
        }

    for(int i=0;i<front.size();i++)
    {
        double temp=0.0;
        for(int j=0;j<distance.length;j++)
        {
            if(i!=j){
            temp+=Math.pow(distance[i][j]-(dmaxx-dminn),2);
            }
        }
        temp/=distance.length-1;
        temp=Math.sqrt(temp);
        temp*=-1;
        temp+=(double) front.get(i).getCrowdingDistance();
        front.get(i).setCrowdingDistance(temp);
    }
    
    //int k = numberOfObjectives la solucion 0 es ella misma
    for (int i = 0; i < distance.length; i++) {
      Arrays.sort(distance[i]);
      double kDistance=0.0;
      for(int k=numberOfObjectives;k>0;k--)
      {
          kDistance += (dmaxx-dminn) / (distance[i][k]+((dminn)));
      }
      double temp=(double) front.get(i).getCrowdingDistance();
      temp-=kDistance;
      front.get(i).setCrowdingDistance(temp);
    }
  
    
  } // crowdingDistanceAssing            

  public double [][] distanceMatrix_normalize(SolutionSet solutionSet, double maxs[], double mins[]) {
    Solution solutionI, solutionJ;

    //The matrix of distances
    double [][] distance = new double [solutionSet.size()][solutionSet.size()];
    //-> Calculate the distances
    for (int i = 0; i < solutionSet.size(); i++){
      distance[i][i] = 0.0;
      solutionI = solutionSet.get(i);
      for (int j = i + 1; j < solutionSet.size(); j++){
        solutionJ = solutionSet.get(j);
        distance[i][j] = this.distanceBetweenObjectives_normalize(solutionI,solutionJ,maxs,mins);
        distance[j][i] = distance[i][j];
      } // for
    } // for        

    //->Return the matrix of distances
    return distance;
  } // distanceMatrix

  public double distanceBetweenObjectives_normalize(Solution solutionI, Solution solutionJ,double maxs[],double mins[]){
    double diff;    //Auxiliar var
    double distance = 0.0;
    //-> Calculate the euclidean distance
    for (int nObj = 0; nObj < solutionI.getNumberOfObjectives();nObj++){
      diff = (solutionI.getObjective(nObj)/maxs[nObj]-mins[nObj]) - (solutionJ.getObjective(nObj)/maxs[nObj]-mins[nObj]);
      distance += Math.pow(diff,2.0);
    } // for  
      return distance;
    }
  


} // Distance

