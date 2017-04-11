package ie.errity.pd.genetic;

import ie.errity.pd.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.BitSet;
import java.util.Random;
import java.lang.Math;
import java.util.Arrays;


/**
 * Provides a means of evolving a population of 
 *{@link  ie.errity.pd.Prisoner Prisoners} 
 *via a genetic algorithm
 * @author	Andrew Errity 99086921
 * @author      Sherri Goings (modified 4/4/2017)
 */
public class Breeder extends JPanel
{
    private Prisoner curPopulation[];
    private double mutateP, crossP; //mutate, cross probabilities
    private int selection, selParam, seed; // selection method, any associated selection parameter, random number seed
    private int popSize;
    private Random rand;
	
    /**
     *Create a new Genetic Breeder
     */
    public Breeder()
    {	
	//defaults
	mutateP = .001;
	crossP = .95;
	selection = 0;
	selParam = 1;
	seed = -1;
	rand = new Random(); //uses current time as seed
    }
	
    /**
     *Create a new Genetic Breeder using {@link  ie.errity.pd.Rules Rules} given
     *@param r1 parameters for the breeder
     */
    public Breeder(Rules r1)
    {	
	//init based on rules	
	mutateP = r1.getMutateP();
	crossP = r1.getCrossP();
	selection = r1.getSelection();
	selParam = r1.getSelectionParameter();
	seed = r1.getSeed();
	if (seed==-1)
	    rand = new Random();
	else
	    rand = new Random(seed); //uses seed set in rules for repeatability
    }

    private int getScoreSum(Prisoner prisoners[]) {
        int scoreSum = 0;
        for (Prisoner p: prisoners) {
            scoreSum += p.getScore();
        }
        return scoreSum;
    }

    private double[] getScaledFitness(Prisoner prisoners[], int scoreSum) {
        double scaledFitnesses[] = new double[popSize];
        double mean = scoreSum / popSize;

        double sumSquaredErr = 0;
        for (Prisoner p: prisoners) {
            sumSquaredErr += Math.pow(p.getScore() - mean, 2);
        }

        double stDev = Math.sqrt(sumSquaredErr / popSize);
        if (stDev == 0) {
            for (int i = 0; i < popSize; i++) {
                double scaledFitness = 1 + (prisoners[i].getScore() - mean) / (2 * stDev);
                scaledFitnesses[i] = scaledFitness < .1 ? .1 : scaledFitness;
            }
        } else {
            for (int i = 0; i < popSize; i++) {
                scaledFitnesses[i] = 1;
            }
        }
        return scaledFitnesses;
    }

    private double sum(double arr[])  {
        double sumSofar = 0;
        for (double x : arr)  {
            sumSofar += x;
        }
        return sumSofar;
    }

    private void fitnessProportionalSelection(Prisoner[] sampleFrom, Prisoner[] sampleTo) {
        int scoreSum = getScoreSum(sampleFrom);

        double scaledFitnesses[] = getScaledFitness(sampleFrom, scoreSum);
        System.out.println("popSize = " + popSize);
        System.out.println("selParam = " + selParam);
        System.err.println("scaledFitnesses =" + Arrays.toString(scaledFitnesses));

        double intervalSize = sum(scaledFitnesses) / (popSize - selParam);

        double distToNextTick = rand.nextDouble() * intervalSize;
        //TEMP

        int curIndex = 0;
        int ticksLeftToPlace = popSize - selParam;
        int i = 0;

//        System.err.println(Arrays.toString(scaledFitnesses));
        System.err.println("intervalSize: " + intervalSize);
//        System.err.println("starting sampling: " + (popSize - numElite));
        while (ticksLeftToPlace > 0) {
            if (distToNextTick < scaledFitnesses[curIndex]) {
                Prisoner sampledPrisoner = (Prisoner) sampleFrom[curIndex].clone();
                sampleTo[ticksLeftToPlace + selParam - 1] = sampledPrisoner;
                System.err.println("curIndex =" + curIndex);
                System.err.println("counter =" + ticksLeftToPlace);
                System.err.println("scaledFitnesses =" + Arrays.toString(scaledFitnesses));
                System.err.println("sum = " + sum(scaledFitnesses));
                System.err.println("distToNextTick = " + distToNextTick);
                ticksLeftToPlace--;

                scaledFitnesses[curIndex] -= distToNextTick;
                distToNextTick = intervalSize;
            } else {
                distToNextTick -= scaledFitnesses[curIndex];
                curIndex = (curIndex + 1) % popSize;
            }
        }
        return;
    }




    /**
     *Breeds the next generation (panmictic mating) of an array of 
     *{@link  ie.errity.pd.Prisoner Prisoners} 
     *@param c	initial population (raw fitness of population must be calcualted previously)
     *@return the next generation
     */
    public Prisoner[] Breed(Prisoner[] c)
    {
	curPopulation = c;	//population to breed
	popSize = curPopulation.length;
	Prisoner Selected[] = new Prisoner[popSize]; // parent pop after selection 
	
		
	// Select parents for next gen 
	// ***ADD CODE (make separate functions) to perform each of the types of selection***
	// selection is an int that determines the method to use
	// selParam is an int that can be used as a parameter for a selection method if required
	
	// One silly selection method which uses rand and the selParam is given below
	// Be sure to use "rand" class variable here as your random number generator (i.e. don't create a new one!)

	// selection method 0, use parameter as a threshhold percentage.  If a random number is less than param, 
	// select the best individual from the population.  Otherwise select a random individual.
	// In this method a param of 0 gives random selection, and a param of 100 gives best-wins-all selection.
	switch (selection) {
        // Initial Threshold Random Thingy Sherri Made
	    case 0:
            // find index of most fit individual
            double maxFit = 0;
            int indexBest = 0;
            for (int i=0; i<popSize; i++) {
            if (curPopulation[i].getScore() > maxFit) {
                maxFit = curPopulation[i].getScore();
                indexBest = i;
            }
            }

            // select according to description above for this method
            for (int i=0; i<popSize; i++) {
            int selIndex = 0;
            if (rand.nextInt(100) < selParam)  // rand less than param, select best individual
                {
                selIndex = indexBest;
                }
            else  // otherwise select random individual
                {
                selIndex = rand.nextInt(popSize);
                }
            Selected[i] = (Prisoner)curPopulation[selIndex].clone();
            }
	        break;

        // Fitness Proportional Selection with sigma scaling and stochastic universal sampling, plus optional elitism
        case 1:
            if (selParam > 0) {
                Prisoner best = c[0];
                for (Prisoner p: c) {
                    if (p.getScore() > best.getScore()) {
                        best = p;
                    }
                }
                for (int i = 0; i < selParam; i++) {
                    Selected[i] = (Prisoner) best.clone();
                }
            }
            fitnessProportionalSelection(c, Selected);



            break;
        // tournament selection
        case 2:




            break;

        // otherwise...
        default: // any other selection method fill pop with always cooperate
            for (int i=0; i<popSize; i++)
            Selected[i] = new Prisoner("ALLC");
            break;
	}

			
	//Crossover & Mutate each pair of selected parents	
	BitSet Offspring[] = new BitSet[2];  // temporarily holds 2 children during crossover/mutation
	for (int d=0; d<popSize; d+=2) {
	    // in case of odd population, just mutate and replace last individual
	    if (d+1 >= popSize) {
		Offspring[0] = Genetic.mutate(Selected[d].getStrat(), mutateP, rand);
		Selected[d] = new Prisoner(Offspring[0]);
	    }
	    else {
							
		if(rand.nextDouble() <= crossP) //Cross Over
		    Offspring = Genetic.crossover(Selected[d].getStrat(),Selected[d+1].getStrat(), rand);
		else //clones
		    {
			Offspring[0] = (BitSet)Selected[d].getStrat().clone();
			Offspring[1] = (BitSet)Selected[d+1].getStrat().clone();
		    }
			
		//Mutation
		Offspring[0] = Genetic.mutate(Offspring[0],mutateP, rand);
		Offspring[1] = Genetic.mutate(Offspring[1],mutateP, rand);
			
		//Replacement - we are done with parents d & d+1, so just replace with children without
		// creating an entire new array
		Selected[d] = new Prisoner(Offspring[0]);
		Selected[d+1] = new Prisoner(Offspring[1]);
	    }
	}
	// pass on children pop to be parents of next gen
	curPopulation = Selected;
	repaint();	//update display (if any)
	return curPopulation; //return the bred population
    }
	
	
    /**
     *Responsible for updating the graphical representation
     */
    public void paintComponent(Graphics g) 
    {
        super.paintComponent(g); //paint background

	//Get limits of viewable area
      	Insets insets = getInsets();
      	int x0 = insets.left;
      	int y0 = insets.top;
	int currentWidth = getWidth() - insets.left - insets.right;
	int currentHeight = getHeight() - insets.top - insets.bottom;
	    
	//Display a series of rectangles, representing the players
	for(int i = 0; i < popSize; i++)
	    {
	    	g.setColor(curPopulation[i].getColor());	
	    	g.fillRect((x0*2)+((currentWidth/popSize)*(i)),(currentHeight/4)+y0,(currentWidth/popSize),currentHeight/2);
	    }
    }
	
    /**
     *Set the {@link  ie.errity.pd.Rules Rules}
     *@param r new {@link  ie.errity.pd.Rules Rules}
     */
    public void setRules(Rules r)
    {
	mutateP = r.getMutateP();
	crossP = r.getCrossP();
	selection = r.getSelection();
	selParam = r.getSelectionParameter();
	seed = r.getSeed();
	if (seed > -1)
	    rand.setSeed(seed);
    }
	
    /**
     *Reset the breeder
     */
    public void clear()
    {
	popSize = 0;
	repaint();
    }
}
	
