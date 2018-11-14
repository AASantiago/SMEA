//  NSGAII_main.java
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

package jmetal.metaheuristics.smea;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.ProblemFactory;
import jmetal.problems.ZDT.ZDT3;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
import jmetal.problems.GLT.GLT1;
import jmetal.problems.GLT.GLT2;
import jmetal.problems.GLT.GLT3;
import jmetal.problems.GLT.GLT4;
import jmetal.problems.GLT.GLT5;
import jmetal.problems.GLT.GLT6;
import jmetal.problems.LZ09.LZ09_F1;
import jmetal.problems.LZ09.LZ09_F2;
import jmetal.problems.LZ09.LZ09_F3;
import jmetal.problems.LZ09.LZ09_F4;
import jmetal.problems.LZ09.LZ09_F5;
import jmetal.problems.LZ09.LZ09_F6;
import jmetal.problems.LZ09.LZ09_F7;
import jmetal.problems.LZ09.LZ09_F8;
import jmetal.problems.LZ09.LZ09_F9;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG5;
import jmetal.problems.WFG.WFG6;
import jmetal.problems.WFG.WFG7;
import jmetal.problems.WFG.WFG8;
import jmetal.problems.WFG.WFG9;
import jmetal.problems.ZDT.ZDT1;
import jmetal.problems.ZDT.ZDT2;
import jmetal.problems.ZDT.ZDT4;
import jmetal.problems.ZDT.ZDT5;
import jmetal.problems.ZDT.ZDT6;
import jmetal.problems.cec2009Competition.UF1;
import jmetal.problems.cec2009Competition.UF10;
import jmetal.problems.cec2009Competition.UF2;
import jmetal.problems.cec2009Competition.UF3;
import jmetal.problems.cec2009Competition.UF4;
import jmetal.problems.cec2009Competition.UF5;
import jmetal.problems.cec2009Competition.UF6;
import jmetal.problems.cec2009Competition.UF7;
import jmetal.problems.cec2009Competition.UF8;
import jmetal.problems.cec2009Competition.UF9;

/** 
 * Class to configure and execute the NSGA-II algorithm.  
 *     
 * Besides the classic NSGA-II, a steady-state version (ssNSGAII) is also
 * included (See: J.J. Durillo, A.J. Nebro, F. Luna and E. Alba 
 *                  "On the Effect of the Steady-State Selection Scheme in 
 *                  Multi-Objective Genetic Algorithms"
 *                  5th International Conference, EMO 2009, pp: 183-197. 
 *                  April 2009)
 */ 

public class SMEA_main {
  public static Logger      logger_ ;      // Logger object
  public static FileHandler fileHandler_ ; // FileHandler object

  /**
   * @param args Command line arguments.
   * @throws JMException 
   * @throws IOException 
   * @throws SecurityException 
   * Usage: three options
   *      - jmetal.metaheuristics.nsgaII.NSGAII_main
   *      - jmetal.metaheuristics.nsgaII.NSGAII_main problemName
   *      - jmetal.metaheuristics.nsgaII.NSGAII_main problemName paretoFrontFile
   */
  public static void main(String [] args) throws 
                                  JMException, 
                                  SecurityException, 
                                  IOException, 
                                  ClassNotFoundException {
    long start = System.currentTimeMillis();
    Problem   problem = null   ; // The problem to solve
    Algorithm algorithm ; // The algorithm to use
    Operator  selection ; // Selection operator
    
    HashMap  parameters ; // Operator parameters
    
    QualityIndicator indicators ; // Object to get quality indicators

    // Logger object and file to store log messages
    logger_      = Configuration.logger_ ;
    fileHandler_ = new FileHandler("SMEA_main.log"); 
    logger_.addHandler(fileHandler_) ;
        
    int poblacion=200;
    int evaluaciones=45000;
    indicators = null ;
    if(args.length==0){
      // Default problem
      //problem = new UF1("Real");  
      //problem = new Kursawe("Real", 3); 
      //problem = new Kursawe("BinaryReal", 3);
      //problem = new Water("Real");
      //problem = new LZ09_F4("Real",21,1,24);
      //problem = new LZ09_F5("Real",21,1,26);
      //problem = new LZ09_F9("Real",21,1,22);
      //problem =new LZ09_F7("Real",21,3,21);
      problem = new GLT5("Real");
      //problem = new ZDT1("Real");
      //problem = new LZ09_F4("Real");
      //problem = new ConstrEx("Real");
      //problem = new DTLZ7("Real");
      //problem = new OKA2("Real") ;
    } else if (args.length == 1) {
      Object [] params = {"Real"};
      problem = (new ProblemFactory()).getProblem(args[0],params);
    } // if
    else if (args.length == 2) {
      Object [] params = {"Real"};
      problem = (new ProblemFactory()).getProblem(args[0],params);
      indicators = new QualityIndicator(problem, args[1]) ;
    } // if
    else if(args[0].compareTo("WFG")==0) {//WFG
        System.out.println("WFG");
        int problema=Integer.valueOf(args[1]);
        int valork=Integer.valueOf(args[2]);
        int valorl=Integer.valueOf(args[3]);
        int objetivos=Integer.valueOf(args[4]);
        poblacion=Integer.valueOf(args[5]);
        evaluaciones=Integer.valueOf(args[6]);
        switch(problema)
        {
            case 1:
                problem=new WFG1("Real",valork,valorl,objetivos);
                break;
            case 2:
                problem=new WFG2("Real",valork,valorl,objetivos);
                break;
            case 3:
                problem=new WFG3("Real",valork,valorl,objetivos);
                break;
            case 4:
                problem=new WFG4("Real",valork,valorl,objetivos);
                break;
            case 5:
                problem=new WFG5("Real",valork,valorl,objetivos);
                break;
            case 6:
                problem=new WFG6("Real",valork,valorl,objetivos);
                break;
            case 7:
                problem=new WFG7("Real",valork,valorl,objetivos);
                break;
            case 8:
                problem=new WFG8("Real",valork,valorl,objetivos);
                break;
            case 9:
                problem=new WFG9("Real",valork,valorl,objetivos);
                break;
            default:
                System.out.println("Problema no seleccionado WFG family");
                break;
                
                
        }
    }else if(args[0].compareTo("DTLZ")==0){ //DTLZ
        System.out.println("DTLZ");
        int problema=Integer.valueOf(args[1]);
        int variables=Integer.valueOf(args[2]);
        int objetivos=Integer.valueOf(args[3]);
        poblacion=Integer.valueOf(args[4]);
        evaluaciones=Integer.valueOf(args[5]);
        switch(problema)
        {
            case 1:
                problem=new DTLZ1("Real",variables,objetivos);
                break;
            case 2:
                problem=new DTLZ2("Real",variables,objetivos);
                break;
            case 3:
                problem=new DTLZ3("Real",variables,objetivos);
                break;
            case 4:
                problem=new DTLZ4("Real",variables,objetivos);
                break;
            case 5:
                problem=new DTLZ5("Real",variables,objetivos);
                break;
            case 6:
                problem=new DTLZ6("Real",variables,objetivos);
                break;
            case 7:
                problem=new DTLZ7("Real",variables,objetivos);
                break;
            default:
                System.out.println("Problema no seleccionado DTLZ family");
                break;
                
        }
    
    } else if(args[0].compareTo("LZ09")==0){ //LZ09
        System.out.println("LZ09");
        int problema=Integer.valueOf(args[1]);
        poblacion=Integer.valueOf(args[2]);
        evaluaciones=Integer.valueOf(args[3]);
        switch(problema)
        {
            case 1:
                problem = new LZ09_F1("Real",21,1,21);
                break;
            case 2:
                problem = new LZ09_F2("Real",21,1,22);
                break;
            case 3:
                problem = new LZ09_F3("Real",21,1,23);
                break;
            case 4:
                problem = new LZ09_F4("Real",21,1,24);
                break;
            case 5:
                problem = new LZ09_F5("Real",21,1,26);
                break;
            case 6:
                problem = new LZ09_F6("Real",31,1,32);
                break;
            case 7:
                problem =new LZ09_F7("Real",21,3,21);
                break;
            case 8:
                problem =new LZ09_F8("Real",21,4,21);
                break;
            case 9:
                problem = new LZ09_F9("Real",22,1,22);
                break;
            default:
                System.out.println("Problema no seleccionado ZDT family");
                break;
        }
     
    }else if(args[0].compareTo("ZDT")==0){ //ZDT
        System.out.println("ZDT");
        int problema=Integer.valueOf(args[1]);
        poblacion=Integer.valueOf(args[2]);
        evaluaciones=Integer.valueOf(args[3]);
        switch(problema)
        {
            case 1:
                problem = new ZDT1("Real",30);
                break;
            case 2:
                problem = new ZDT2("Real",30);
                break;
            case 3:
                problem = new ZDT3("Real",30);
                break;
            case 4:
                problem = new ZDT4("Real",10);
                break;
            case 5:
                problem = new ZDT5("Binary",11);
                break;
            case 6:
                problem = new ZDT6("Real",10);
                break;
            default:
                System.out.println("Problema no seleccionado ZDT family");
                break;
        }
     
    }else if(args[0].compareTo("UF")==0){ //LZ09
        System.out.println("UF");
        int problema=Integer.valueOf(args[1]);
        poblacion=Integer.valueOf(args[2]);
        evaluaciones=Integer.valueOf(args[3]);
        switch(problema)
        {
            case 1:
                problem = new UF1("Real");
                break;
            case 2:
                problem = new UF2("Real");
                break;
            case 3:
                problem = new UF3("Real");
                break;
            case 4:
                problem = new UF4("Real");
                break;
            case 5:
                problem = new UF5("Real");
                break;
            case 6:
                problem = new UF6("Real");
                break;
            case 7:
                problem =new UF7("Real");
                break;
            case 8:
                problem =new UF8("Real");
                break;
            case 9:
                problem = new UF9("Real");
                break;
            case 10:
                problem = new UF10("Real");
                break;      
            default:
                System.out.println("Problema no seleccionado CEC09 family");
                break;
        }
    }else if(args[0].compareTo("GLT")==0){
     System.out.println("ZDT");
        int problema=Integer.valueOf(args[1]);
        poblacion=Integer.valueOf(args[2]);
        evaluaciones=Integer.valueOf(args[3]);
        switch(problema)
        {
            case 1:
                problem = new GLT1("Real");
                break;
            case 2:
                problem = new GLT2("Real");
                break;
            case 3:
                problem = new GLT3("Real");
                break;
            case 4:
                problem = new GLT4("Real");
                break;
            case 5:
                problem = new GLT5("Real");
                break;
            case 6:
                problem = new GLT6("Real");
                break;
            default:
                System.out.println("Problema no seleccionado GLT family");
                break;
        }
    }else{
      // Default problem
      System.out.println("No se selecciono ningun problema a resolver");
      //problem = new Kursawe("Real", 3); 
      //problem = new Kursawe("BinaryReal", 3);
      //problem = new Water("Real");
      //problem = new LZ09_F1("Real",21,1,21);
      //problem = new ConstrEx("Real");
      //problem = new DTLZ7("Real");
      //problem = new OKA2("Real") ;
  } // else
    
    algorithm = new SMEA(problem);
    //algorithm = new ssNSGAII(problem);

    // Algorithm parameters
    algorithm.setInputParameter("populationSize",poblacion);
    algorithm.setInputParameter("maxEvaluations",evaluaciones);

    // Selection Operator 
    parameters = null ;
    selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters) ;                           

    // Add the operators to the algorithm
    algorithm.addOperator("selection",selection);

    // Add the indicator object to the algorithm
    algorithm.setInputParameter("indicators", indicators) ;
    
    // Execute the Algorithm
    long initTime = System.currentTimeMillis();
    SolutionSet population = algorithm.execute();
    long estimatedTime = System.currentTimeMillis() - initTime;
    
    // Result messages 
    logger_.info("Total execution time: "+estimatedTime + "ms");
    logger_.info("Variables values have been writen to file VAR");
    population.printVariablesToFile("VAR");    
    logger_.info("Objectives values have been writen to file FUN");
    population.printObjectivesToFile("FUN");
    long end = System.currentTimeMillis();
    long res = end - start;
    System.out.println("Minutos: "+((res/1000)/60));
    if (indicators != null) {
      logger_.info("Quality indicators") ;
      logger_.info("Hypervolume: " + indicators.getHypervolume(population)) ;
      logger_.info("GD         : " + indicators.getGD(population)) ;
      logger_.info("IGD        : " + indicators.getIGD(population)) ;
      logger_.info("Spread     : " + indicators.getSpread(population)) ;
      logger_.info("Epsilon    : " + indicators.getEpsilon(population)) ;  
     
      int evaluations = ((Integer)algorithm.getOutputParameter("evaluations")).intValue();
      logger_.info("Speed      : " + evaluations + " evaluations") ;
    } // if
  } //main
} // NSGAII_main
