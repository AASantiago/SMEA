//  HUXCrossover.java
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

package jmetal.operators.crossover;

import jmetal.core.Solution;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.BinarySolutionType;
import jmetal.encodings.variable.Binary;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import jmetal.encodings.solutionType.ArrayRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.wrapper.XReal;

/**
 * This class allows to apply a HUX crossover operator using two parent
 * solutions.
 * NOTE: the operator is applied to the first encodings.variable of the solutions, and
 * the type of the solutions must be Binary or BinaryReal
 */
public class GeneticDrift extends Crossover{

	/**
   * Valid solution types to apply this operator 
   */
  private static final List VALID_TYPES = Arrays.asList(RealSolutionType.class,
  		                                            ArrayRealSolutionType.class, BinarySolutionType.class) ;

  private Double probability_ = null ;
  /**
   * Constructor
   * Create a new instance of the HUX crossover operator.
   */
  public GeneticDrift(HashMap<String, Object> parameters) {
    super(parameters) ;
    
  	if (parameters.get("probability") != null)
  		probability_ = (Double) parameters.get("probability") ;  		
  } // HUXCrossover


   /**
   * Constructor
   * Create a new intance of the HUX crossover operator.
   */
   //public HUXCrossover(Properties properties) {
   // this();
   //} // HUXCrossover



  /**
   * Perform the crossover operation
   * @param probability Crossover probability
   * @param parent1 The first parent
   * @param parent2 The second parent
   * @return An array containing the two offsprings
   * @throws JMException 
   */
  public Solution[] doCrossover(double   probability, 
                                Solution parent1, 
                                Solution parent2) throws JMException {
    Solution [] offSpring = new Solution[2];
    offSpring[0] = new Solution(parent1);
    offSpring[1] = new Solution(parent2);
    XReal childuno=new XReal(offSpring[0]);
    XReal childdos=new XReal(offSpring[1]);
    try {         
      if (PseudoRandom.randDouble() < probability) {
        for (int var = 0; var < parent1.getDecisionVariables().length; var++) {
          XReal p1 = new XReal(parent1);
          XReal p2 = new XReal(parent2);


          for (int bit = 0; bit < p1.getNumberOfDecisionVariables(); bit++) {
              if (PseudoRandom.randDouble() < 0.5) {
                childuno.setValue(bit,p2.getValue(bit));
                childdos.setValue(bit,p1.getValue(bit));
              }else{
                  childuno.setValue(bit,p1.getValue(bit));
                  childdos.setValue(bit,p2.getValue(bit));
              }
          }
        }  
        
      }          
    }catch (ClassCastException e1) {
      
      Configuration.logger_.severe("HUXCrossover.doCrossover: Cannot perfom " +
          "SinglePointCrossover ") ;
      Class cls = java.lang.String.class;
      String name = cls.getName(); 
      throw new JMException("Exception in " + name + ".doCrossover()") ;
    }        

    return offSpring;                                                                                      
  } // doCrossover

  
  /**
  * Executes the operation
  * @param object An object containing an array of two solutions 
  * @return An object containing the offSprings
  */
  public Object execute(Object object) throws JMException {
    Solution [] parents = (Solution [])object;
    
    if (parents.length < 2)
    {
      Configuration.logger_.severe("HUXCrossover.execute: operator needs two " +
          "parents");
      Class cls = java.lang.String.class;
      String name = cls.getName(); 
      throw new JMException("Exception in " + name + ".execute()") ;      
    }

    if (!(VALID_TYPES.contains(parents[0].getType().getClass())  &&
        VALID_TYPES.contains(parents[1].getType().getClass())) ) {

      Configuration.logger_.severe("HUXCrossover.execute: the solutions " +
          "are not of the right type. The type should be 'Binary' of " +
          "'BinaryReal', but " +
          parents[0].getType() + " and " + 
          parents[1].getType() + " are obtained");

      Class cls = java.lang.String.class;
      String name = cls.getName(); 
      throw new JMException("Exception in " + name + ".execute()") ;

    } // if 
        
    Solution [] offSpring = doCrossover(probability_,
                                        parents[0],
                                        parents[1]);
    
    for (int i = 0; i < offSpring.length; i++)
    {
      offSpring[i].setCrowdingDistance(0.0);
      offSpring[i].setRank(0);
    } 
    
    return offSpring;
  } // execute
} // HUXCrossover