//  SMEA.java
//
//  Author:
//       Alejandro Santiago <aurelio.santiago@upalt.edu.mx>
//
//  Copyright (c) 2018 Alejandro Santiago
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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import jmetal.core.*;
import jmetal.metaheuristics.moead.Utils;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.Hypervolume;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.Ranking;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.comparators.FitnessComparator;

/** 
 *  Implementation of SMEA.
 */

public class SMEA extends Algorithm {
  /**
   * Constructor
   * @param problem Problem to solve
   */
    
  public double grid[][][]=null;
  public Solution gridsoluciones[][]=null;
  SolutionSet population = null;
  public jmetal.qualityIndicator.util.MetricsUtil utils_;
  public SMEA(Problem problem) {
    super (problem) ;
  } // SMEA


  public SolutionSet execute() throws JMException, ClassNotFoundException {
    int populationSize;
    int maxEvaluations;
    int evaluations;

    QualityIndicator indicators; // QualityIndicator object
    int requiredEvaluations; // Use in the example of use of the
    // indicators object (see below)

    SolutionSet offspringPopulation;
    SolutionSet union;

    Operator mutationOperator;
    Operator crossoverOperator;
    Operator selectionOperator;
    Operator crossoverOperator_DE;
    HashMap parameters = new HashMap() ;
    parameters.put("CR", 1.0) ;
    parameters.put("F", 0.9) ;
    crossoverOperator_DE = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);                   
    Operator mutationPolynomial;
    parameters = new HashMap() ;
    parameters.put("probability", 1.0/problem_.getNumberOfVariables()) ;
    parameters.put("distributionIndex", 20.0) ;
    mutationPolynomial = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

    //Read the parameters
    populationSize = ((Integer) getInputParameter("populationSize")).intValue();
    maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();
    indicators = (QualityIndicator) getInputParameter("indicators");

    //Initialize the variables
    //population = new SolutionSet(populationSize);
    SolutionSet nuevapoblacion = null;
    evaluations = 0;

    requiredEvaluations = 0;

    //Read the operators
    mutationOperator = operators_.get("mutation");
    crossoverOperator = operators_.get("crossover");
    selectionOperator = operators_.get("selection");
    Solution newSolution;
    
    if(problem_.getNumberOfObjectives()==2)
    {
        populationSize=1*100;
        grid=new double[1][100][problem_.getNumberOfVariables()];
        gridsoluciones=new Solution[1][100];
    }
    if(problem_.getNumberOfObjectives()==3)
    {
        populationSize=7*15;
        grid=new double[7][15][problem_.getNumberOfVariables()];
        gridsoluciones=new Solution[7][15];
        
    }
    grid=new double[20][10][problem_.getNumberOfVariables()];
    gridsoluciones=new Solution[20][10];
    
    double learning_rate=0.7;
    double tamanio_vecindario=0.5*Math.sqrt((Math.pow(grid.length,2)+Math.pow(grid[0].length,2))/populationSize); //corregido tipo zhang
    // Generations
    int generations=0;
    int maxgenerations=maxEvaluations/populationSize; //corregido tipo zhang
    
    double[][] distMatrix=new double[populationSize][populationSize];
    for(int x=0;x<populationSize;x++)
    {
        for(int j=0;j<populationSize;j++)
        {
            distMatrix[x][j]=Math.sqrt(Math.pow(x-j,2)); //corregido tipo zhang (matriz simetrica)
        }
    }
    
    int neigDepot[][]=new int[populationSize][5];//busca los 5 indices mas cercanos
    for(int x=0;x<populationSize;x++)
    {
        int indices[]=new int[populationSize];
        double distancias[]=new double[populationSize];
        for(int j=0;j<populationSize;j++)
        {
            indices[j]=j;
            distancias[j]=distMatrix[x][j];//matriz simetrica
        }
        //SORT!!!
        quicksort(distancias,indices,0,populationSize-1);
        for(int j=0;j<5;j++)
        neigDepot[x][j]=indices[j+1];
    }
    /*
    for(int x=0;x<populationSize;x++)
    {
        for(int j=0;j<5;j++)
        {
            System.out.print(neigDepot[x][j]+"\t");
        }
        System.out.println();
    }*/
    population=new SolutionSet(populationSize);
    for(int x=0;x<populationSize;x++)
    {
        newSolution = new Solution(problem_);
        problem_.evaluate(newSolution);
        problem_.evaluateConstraints(newSolution);
        evaluations++;
        population.add(newSolution);
    }
    

      
      SolutionSet traindata=new SolutionSet(populationSize);
      SolutionSet archive=new SolutionSet(populationSize);
      for(int x=0;x<populationSize;x++)
      {
          traindata.add(population.get(x));
          archive.add(population.get(x));
      }
     int[] permutation = new int[populationSize];
     Utils.randomPermutation(permutation,populationSize);
 
      double pesos[][]=new double[populationSize][problem_.getNumberOfVariables()];
      //distribuye los pesos originales aleatoriamente
      for(int j=0;j<populationSize;j++)
      for(int x=0;x<problem_.getNumberOfVariables();x++) //corregido tipo zhang
          pesos[j][x]=population.get(permutation[j]).getDecisionVariables()[x].getValue();
      
      double tau0=0.7;
      double rad0=Math.sqrt((Math.pow(grid.length,2)+Math.pow(grid[0].length , 2)))/2;
      for(generations=1;generations<=maxgenerations;generations++)
      {
        System.out.println("SMEA Generacion "+generations);
       if(traindata.size()>0)//si hay patrones para entrenar
       {
           
            SOMOperator(tau0,rad0,traindata,generations,maxEvaluations,populationSize,distMatrix,pesos);
            
           
       }
       for(int x=0;x<populationSize;x++)
       {
           Solution[] offSpring=new Solution[2];
           Solution parents[]=new Solution[5];
           if(Math.random()<0.9)
           {
               BinaryTournament(neigDepot[x],parents);
               parents[2]=population.get(x);
               parents[3]=population.get(x);
               parents[4]=population.get(x);
           }else{
               permutation = new int[populationSize];
               Utils.randomPermutation(permutation,populationSize);
               BinaryTournament(permutation,parents);
               parents[2]=population.get(x);
               parents[3]=population.get(x);
               parents[4]=population.get(x);
           }
           
           //DE
           offSpring[0]=(Solution) crossoverOperator_DE.execute(new Object[]{population.get(x), parents});
           mutationPolynomial.execute(offSpring[0]);
           problem_.evaluate(offSpring[0]);
           problem_.evaluateConstraints(offSpring[0]);
           Removesolution(offSpring[0],archive);
       }
          
       traindata.clear();
       int flag=-1;
       for(int x=0;x<archive.size();x++)
       {
           flag=1;
           for(int j=0;j<population.size();j++)
           {
               for(int z=0;z<problem_.getNumberOfVariables();z++)
               {
                   
                   if(population.get(x).getDecisionVariables()[z].getValue()!=archive.get(j).getDecisionVariables()[z].getValue())
                   {
                       flag=0;
                   }
               }
           }
           if(flag==1)
               traindata.add(archive.get(x));
               
       }
       population=archive;
          
      }

    // Return as output parameter the required evaluations
    setOutputParameter("evaluations", requiredEvaluations);

    // Return the first non-dominated front
    Ranking ranking = new Ranking(population);
    ranking.getSubfront(0).printFeasibleFUN("FUN_NSGAII") ;

    return ranking.getSubfront(0);
  } // execute
  
  public void Removesolution(Solution nueva, SolutionSet soluciones)
  {
      int remain=soluciones.size();
      int index=0;
      soluciones.setCapacity(soluciones.size()+1);
      soluciones.add(nueva);
      Ranking ranking=new Ranking(soluciones);
      soluciones.clear();
      SolutionSet front = ranking.getSubfront(0);
      Distance distance = new Distance();
      if(front.size()==soluciones.size())
      {
          //POR HYPERVOLUMEN
          Hypervolume ab=new Hypervolume();
          double [] maximumValues = utils_.getMaximumValues(front.writeObjectivesToMatrix(),problem_.getNumberOfObjectives());
          double [] minimumValues = utils_.getMinimumValues(front.writeObjectivesToMatrix(),problem_.getNumberOfObjectives());
          // STEP 2. Get the normalized front
          double [][]normalizedFront = utils_.getNormalizedFront(front.writeObjectivesToMatrix(), maximumValues, minimumValues);
          double[][] invertedFront = utils_.invertedFront(normalizedFront);
          double hv_original=ab.calculateHypervolume(normalizedFront, front.size(),problem_.getNumberOfObjectives());
          double peorcontribucion=Double.MAX_VALUE;
          int peorsol=-1;
          int contee=0;
          for(int x=0;x<front.size();x++){
          double [][]tempfront=new double[front.size()-1][problem_.getNumberOfObjectives()];
          int cont_tempfront=0;
          for(int k=0;k<front.size();k++)
          {
              for(int q=0;q<problem_.getNumberOfObjectives();q++)
              {
                  if(k!=x)
                  tempfront[cont_tempfront][q]=front.get(k).getObjective(q);
              }
              if(k!=x)
              cont_tempfront++;
          }
          // STEP 2. Get the normalized front
          normalizedFront = utils_.getNormalizedFront(tempfront, maximumValues, minimumValues);
          invertedFront = utils_.invertedFront(normalizedFront);
          double nuevo=ab.calculateHypervolume(invertedFront, front.size()-1, problem_.getNumberOfObjectives());
          if(hv_original-nuevo<peorcontribucion)
          {
              peorcontribucion=hv_original-nuevo;
              peorsol=x;
              
          }
          
          }
          front.remove(peorsol);
          soluciones.clear();
          for(int x=0;x<front.size();x++)
          {
              soluciones.add(front.get(x));
          }
      }else
      {
          while ((remain > 0) && (remain >= front.size())) {
        //Assign crowding distance to individuals
        distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
        //Add the individuals of this front
        for (int k = 0; k < front.size(); k++) {
          soluciones.add(front.get(k));
        } // for

        //Decrement remain
        remain = remain - front.size();

        //Obtain the next front
        index++;
        if (remain > 0) {
          front = ranking.getSubfront(index);
        } // if        
      } // while

      // Remain is less than front(index).size, insert only the best one
      if (remain > 0) {  // front contains individuals to insert                        
        distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
        front.sort(new CrowdingComparator());
        for (int k = 0; k < remain; k++) {
          soluciones.add(front.get(k));
        } // for
      }
      }
      
  }
  public void BinaryTournament(int vecindario[],Solution parents[])
  {
      int pos1,pos2;
      pos1=pos2=-1;
      Random rnd=new Random();
      DominanceComparator comparador=new DominanceComparator();
      while(pos1==pos2)
      {
          pos1=rnd.nextInt(vecindario.length);
          pos2=rnd.nextInt(vecindario.length);
      }
      
      if(comparador.compare(population.get(vecindario[pos1]),population.get(vecindario[pos2]))==-1)
      {
          parents[0]=population.get(vecindario[pos1]);
      }else if(comparador.compare(population.get(vecindario[pos2]),population.get(vecindario[pos1]))==-1)
      {
          parents[0]=population.get(vecindario[pos2]);
      }else{
      
          if(Math.random()<0.5)
              parents[0]=population.get(vecindario[pos1]);
          else
              parents[0]=population.get(vecindario[pos2]);
      }
      
      pos1=pos2=-1;
      while(pos1==pos2)
      {
          pos1=rnd.nextInt(vecindario.length);
          pos2=rnd.nextInt(vecindario.length);
      }
      
      if(comparador.compare(population.get(vecindario[pos1]),population.get(vecindario[pos2]))==-1)
      {
          parents[1]=population.get(vecindario[pos1]);
      }else if(comparador.compare(population.get(vecindario[pos2]),population.get(vecindario[pos1]))==-1)
      {
          parents[1]=population.get(vecindario[pos2]);
      }else{
      
          if(Math.random()<0.5)
              parents[1]=population.get(vecindario[pos1]);
          else
              parents[1]=population.get(vecindario[pos2]);
      }
      
      
      
  }
  public void SOMOperator(double tau0, double rad0, SolutionSet traindata, int generations, int maxEvaluations, int populationSize, double[][] distMatrix, double[][] pesos) throws JMException
  {
           int [] permutation = new int[traindata.size()];
           Utils.randomPermutation(permutation, traindata.size());
           for(int x=0;x<traindata.size();x++)//por cada patron a entrenar
           {
            double tau=tau0*(1-((generations-1)*populationSize+x)/maxEvaluations);
            double rad=rad0*(1-((generations-1)*populationSize+x)/maxEvaluations);
            Solution temp=new Solution(traindata.get(permutation[x]));//ESTO NO VENIA ALEATORIO EN EL CODIGO ORIGINAL
            //WINUNIT
            double win_distancia=Double.MAX_VALUE;
            int win_pos=-1;
            for(int j=0;j<populationSize;j++)
            {
                double distancia=0;
                for(int w=0;w<problem_.getNumberOfVariables();w++)
                {
                    distancia+=Math.pow(temp.getDecisionVariables()[w].getValue()-pesos[j][w],2);
                }
                distancia=Math.sqrt(distancia);
                if(distancia<win_distancia)
                {
                    win_distancia=distancia;
                    win_pos=j;
                }
            }
            //WEIGHT UPDATE
            for(int j=0;j<populationSize;j++)
            {
                
                if(distMatrix[win_pos][j]<rad)
                {
                    for(int w=0;w<problem_.getNumberOfVariables();w++)
                    pesos[j][w]+=tau*Math.exp(-1*distMatrix[win_pos][j])*temp.getDecisionVariables()[w].getValue()-pesos[j][w];
                }
            }
              
           }
                                  
            //CLUSTERING PROCESS 
           permutation = new int[population.size()];
           double pesostemp[][]=new double [population.size()][problem_.getNumberOfVariables()];
           for(int j=0;j<populationSize;j++)
           {for(int x=0;x<problem_.getNumberOfVariables();x++)
               pesostemp[j][x]=pesos[j][x];
           }
           Utils.randomPermutation(permutation, populationSize);
           int winunits[]=new int[populationSize];
           for(int x=0;x<populationSize;x++)
           {
                double best_distancia=Double.MAX_VALUE;
                int best_pos=-1;
                for(int j=0;j<populationSize;j++)
                {
                    double distancia=0;
                    for(int w=0;w<problem_.getNumberOfVariables();w++)
                    {
                           distancia+=Math.pow((pesostemp[j][w]-population.get(permutation[x]).getDecisionVariables()[w].getValue()),2);
                    }
                    distancia=Math.sqrt(distancia);
                    if(distancia<best_distancia)
                    {
                        best_distancia=distancia;
                        best_pos=j;
                    }
                }
                winunits[permutation[x]]=best_pos;
                for(int w=0;w<problem_.getNumberOfVariables();w++)
                {
                    pesostemp[best_pos][w]=Double.MAX_VALUE;
                }
           }
           
           //por eso se reacomoda la población en sus indices jejeje..
           SolutionSet original=new SolutionSet(populationSize);
           for(int x=0;x<populationSize;x++)
               original.add(population.get(x));
           
           for(int x=0;x<populationSize;x++)
           {
               population.replace(x,original.get(winunits[x]));
           }
           
           
  }
    
  public static void quicksort(double valores[], int indices[], int izq, int der) {

  double pivote=valores[izq]; // tomamos primer elemento como pivote
  int pivote2=indices[izq];
  int i=izq; // i realiza la búsqueda de izquierda a derecha
  int j=der; // j realiza la búsqueda de derecha a izquierda
  double aux;
  int aux2;
 
  while(i<j){            // mientras no se crucen las búsquedas
     while(valores[i]<=pivote && i<j) i++; // busca elemento mayor que pivote
     while(valores[j]>pivote) j--;         // busca elemento menor que pivote
     if (i<j) {                      // si no se han cruzado                      
         aux= valores[i];                  // los intercambia
         aux2=indices[i];
         valores[i]=valores[j];
         indices[i]=indices[j];
         valores[j]=aux;
         indices[j]=aux2;
     }
   }
   valores[izq]=valores[j]; // se coloca el pivote en su lugar de forma que tendremos
   indices[izq]=indices[j];
   valores[j]=pivote; // los menores a su izquierda y los mayores a su derecha
   indices[j]=pivote2;
   if(izq<j-1)
      quicksort(valores,indices,izq,j-1); // ordenamos subarray izquierdo
   if(j+1 <der)
      quicksort(valores,indices,j+1,der); // ordenamos subarray derecho
}
} // NSGA-II


