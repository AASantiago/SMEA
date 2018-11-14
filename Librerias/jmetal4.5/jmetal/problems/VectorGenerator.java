/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package jmetal.problems;
import jmetal.core.Solution;
import jmetal.core.Problem;
import jmetal.encodings.solutionType.ArrayRealSolutionType;
import jmetal.util.JMException;

/**
 *
 * @author Alx
 */
//LA POBLACION SE LE TIENE QUE MOVER DIRECTAMENTE EN ARRAYREALSOLUTIONTYPE
public class VectorGenerator extends Problem {
    
    public VectorGenerator(String solutionType) throws ClassNotFoundException{
        this(solutionType,2,1.0);
    }
    
    public VectorGenerator(String solutionType, Integer objetivos, Double suma)
    {
        numberOfVariables_  = objetivos;
        numberOfObjectives_ = 2;
        numberOfConstraints_ = 1;
        problemName_        = "VectorGenerator";
        lowerLimit_ = new double[numberOfVariables_];
        upperLimit_ = new double[numberOfVariables_];
        for(int i=0;i<numberOfVariables_;i++)
        {
            lowerLimit_[i]=0.0;
            upperLimit_[i]=suma;
        }
        
        if(solutionType.compareTo("ArrayReal")==0)
        {
            solutionType_ = new ArrayRealSolutionType(this);
            
        }else{
            System.out.println("Error: solution type " + solutionType + " invalid") ;
            System.exit(-1) ;
        }
        
        
    }

    @Override
    public void evaluate(Solution solution) throws JMException {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public void evaluateConstraints(Solution solution) throws JMException {
    
    }
    
}
