//  BinaryTournament.java
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

package jmetal.operators.selection;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.DominanceComparator;

import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import jmetal.util.comparators.BinaryTournamentComparator;
import jmetal.util.comparators.CrowdingDistanceComparator;
import jmetal.util.comparators.RankAndCrowdingComparator;

/**
 * This class implements an binary tournament selection operator
 */
public class TournamentCrowded extends Selection {

  /**
   * Stores the <code>Comparator</code> used to compare two
   * solutions
   */
  private Comparator comparator_;
  int tournament_size;
  /**
   * Constructor
   * Creates a new Binary tournament operator using a BinaryTournamentComparator
   */
  public TournamentCrowded(HashMap<String, Object> parameters){
  	super(parameters) ;
  	tournament_size = (int) parameters.get("size") ;  	
      comparator_ = new RankAndCrowdingComparator();
  } // BinaryTournament


  
  /**
  * Performs the operation
  * @param object Object representing a SolutionSet
  * @return the selected solution
  */
  public Object execute(Object object){
    //ES SENSIBLE A ESTO EL TORNEO TIENE QUE ESTAR SUPER BIEN HECHO
    //TIPO JMETAL 5
    SolutionSet solutionSet = (SolutionSet)object;
    
    Solution bestKnown,candidate;
    int flag;
    if (solutionSet.size() == 1) {
      bestKnown = solutionSet.get(0);
    } else {
      bestKnown = solutionSet.get(PseudoRandom.randInt(0,solutionSet.size()-1));
      int count = 1; // at least 2 solutions are compared
      do {
       do{ 
       candidate = solutionSet.get(PseudoRandom.randInt(0,solutionSet.size()-1));
       }while(bestKnown.equals(candidate));
        flag = comparator_.compare(bestKnown, candidate);
      if (flag == 1) {
        bestKnown = candidate;
      }
      } while (++count < tournament_size);
    }
    
  return bestKnown;
  } // execute
} // BinaryTournament
