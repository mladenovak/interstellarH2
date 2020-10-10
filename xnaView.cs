using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using System.Windows.Forms;
using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework.Content;
using System.Diagnostics;
using System.ComponentModel;
using System.IO;

namespace InterstellarIce
{
    class xnaView : GraphicsDeviceControl
    {
        ContentManager Content;
		
		public WinFormCamera camera;

		//Model kocka;

		public BackgroundWorker bw;
		public Timer timer;
		//Stopwatch stopwatch;


        /// <summary>
        /// Initializes the control.
        /// </summary>
        protected override void Initialize()
        {
			Content = new ContentManager(Services, "ContentXNA");
			
			camera = new WinFormCamera(
				this.GraphicsDevice,
				 0.7731419f,
				 -0.2597786f,
				0,
				250,
				new Vector3(0f, -10f, 0f));
	
           
			//font = Content.Load<SpriteFont>("SpriteFont1");

			//kocka = Content.Load<Model>("cube");

			instancedModel_grain = Content.Load<Model>("cube_origin");
			instancedModelBones_grain = new Matrix[instancedModel_grain.Bones.Count];
			instancedModel_grain.CopyAbsoluteBoneTransformsTo(instancedModelBones_grain);
			
			instancedModel_ad = Content.Load<Model>("sphere_origin_little");
			instancedModelBones_ad = new Matrix[instancedModel_ad.Bones.Count];
			instancedModel_ad.CopyAbsoluteBoneTransformsTo(instancedModelBones_ad);

			instancedModel_h2 = Content.Load<Model>("sphere2_origin_little");
			instancedModelBones_h2 = new Matrix[instancedModel_h2.Bones.Count];
			instancedModel_ad.CopyAbsoluteBoneTransformsTo(instancedModelBones_h2);


            instanceModelEffect = Content.Load<Effect>("InstancedModelColor");

            Vector3 LightDirection = new Vector3(-4, -10, -5);
            LightDirection.Normalize();
            instanceModelEffect.Parameters["LightDirection"].SetValue(LightDirection);
            instanceModelEffect.Parameters["DiffuseLight"].SetValue(0.7f);
            instanceModelEffect.Parameters["AmbientLight"].SetValue(0.3f);


            grInitializeGrid();
            grCalculateTransitionProbabilities();

            adInitialize();


            // Start the animation timer.
            //stopwatch = Stopwatch.StartNew();

			timer = new Timer();
			timer.Interval = 16;
			timer.Tick += new EventHandler(timer_Tick);
			//timer.Start();

			bw = new BackgroundWorker();
			bw.WorkerSupportsCancellation = true;
			bw.WorkerReportsProgress = true;
			bw.DoWork += new DoWorkEventHandler(bw_DoWork);

            // Hook the idle event to constantly redraw our animation.
            Application.Idle += delegate { Invalidate(); };
        }

		void timer_Tick(object sender, EventArgs e)
		{

			UpdateAll();

		}
		void bw_DoWork(object sender, DoWorkEventArgs e)
		{
			BackgroundWorker worker = sender as BackgroundWorker;
			
			while (!worker.CancellationPending)
			{
				UpdateAll();
			}
		}

		/**** Surface roughening ****/

        private int lx = 1; // linear x dimension
        private int ly = 1; // linear y dimenion

        private float grTemperature = 0.1f; // temperature
        private float grField = 0.1f; // chemical potential difference
        private float grJ = 1f;
        private const int z = 4; // koordinacijski broj za 2D kvadratnu: z=4
        private int[,] grHeigthGrid; // solid heigth
        private float[] grTransitionProbability = new float[6]; // transition probability: k=1 adsorp, k=1-6 desorp
        private int grMag = 0; // magnitude of heigth sum
        private float grNnex = 0; // nearest neighbour energy
		private long grTimeStepCount = 0;


        public int LengthX
        {
            get { return lx; }
            set
            {
                lx = value > 0 ? value : 1;
                grInitializeGrid();
				adInitialize();
            }
        }
        public int LengthY
        {
            get { return ly; }
            set
            {
                ly = value > 0 ? value : 1;
                grInitializeGrid();
				adInitialize();
            }
        }

		public int[,] GrHeigthGrid
		{
			get { return grHeigthGrid; }
		}
        public float GrTemperature
        {
            get { return grTemperature; }
            set
            {
                grTemperature = value > 0 ? value : 0;
                grCalculateTransitionProbabilities();
            }
        }
        public float GrJenergy
        {
            get { return grJ; }
            set
            {
                grJ = value;
                grCalculateTransitionProbabilities();
            }
        }
        public float GrChemicalPotential
        {
            get { return grField; }
            set
            {
				grField = value;
                grCalculateTransitionProbabilities();
            }
        }
        public float[] GrTransitionProbability
        {
            get { return grTransitionProbability; }
        }
		public float GrMeanHeigth
        {
            get { return (float)grMag / (lx * ly); }
        }
        public float GrMeanSurfaceEnergy
        {
            get { return (float)grNnex / (lx * ly); }
        }
		public long GrElapsedTimeSteps
		{
			get { return grTimeStepCount; }
		}

		/*
		public float grGrainMonteCarloFrequency
        {
            get { return mcFrequency; }
            set
            {
                mcFrequency = (value >= 1) ? value : 1;
                mcTimespanGrain = TimeSpan.FromSeconds(1f / mcFrequency);
                elapsedTimeGrain = TimeSpan.Zero;
            }
        }
        public int grGrainMonteCarloLoop
        {
            get { return loop; }
            set { loop = (value >= 1) ? value : 1; }
        }
		*/
        // inicijalizira mrezu
        public void grInitializeGrid()
        {
            grHeigthGrid = new int[lx, ly];

            for (int i = 0; i < lx; i++)
            {
                for (int j = 0; j < ly; j++)
                {
                    //ih[i, j] = (int)MathHelper.Lerp(-2, 2, (float)random.NextDouble());
                    grHeigthGrid[i, j] = 0;
                }
            }
			grTimeStepCount = 0;
			grCalculateEnergy();
			FLAG_grCalculateNextDraw = true;
        }

        // racuna sumu visina i ukupnu energiju susjeda na razlicitim visinama (za J=1)
        public void grCalculateEnergy()
        {
            grMag = 0;
            grNnex = 0;
            int ixRight, iyDown;

            for (int ix = 0; ix < lx; ix++)
            {
                ixRight = (ix < lx - 1) ? ix + 1 : 0; // ciklicki rubni uvjeti

                for (int iy = 0; iy < ly; iy++)
                {
                    iyDown = (iy < ly - 1) ? iy + 1 : 0;

                    grMag += grHeigthGrid[ix, iy];

                    // u svakom kvadratu mreze gleda susjeda ispod sebe i desno od sebe
                    grNnex += grJ * Math.Abs(grHeigthGrid[ix, iy] - grHeigthGrid[ix, iyDown]) + grJ * Math.Abs(grHeigthGrid[ix, iy] - grHeigthGrid[ixRight, iy]);
                }
            }
        }

        // inicijalizira vjerojatnosti prijelaza
        public void grCalculateTransitionProbabilities()
        {
            grTransitionProbability = new float[6];
            // vjerojatnost = w*deltaT <= 1, norma:
			double deltaT = 1 / (Math.Exp(grField / grTemperature) + Math.Exp(z * grJ / grTemperature)); // z=koord broj

            // kristalizacija: ovisi o razlici kemijskog potencijala
			grTransitionProbability[5] = (float)(Math.Exp(grField / grTemperature) * deltaT);

            // evaporacija: vjerojatnosti za razliciti broj susjeda koji su ispod razine
            for (int n = 0; n <= 4; n++)
            {
                double deltaE = -2 * grJ * (n - z / 2);
                grTransitionProbability[n] = (float)(Math.Exp(-deltaE / grTemperature) * deltaT);
            }
            for (int n = 0; n < grTransitionProbability.Length; n++)
            {
                if (float.IsNaN(grTransitionProbability[n]))
                {
                    grTransitionProbability[n] = 0;
                }
            }
        }

        // monte carlo step
        public void grMonteCarloStep()
        {
            // biram nasumicni stupac
            int iX = (int)(Ranlux.NextFloat() * lx);
            if (iX == lx) iX--; // ako random ispadne tono 1.0
            int iY = (int)(Ranlux.NextFloat() * ly);
            if (iY == ly) iY--; // ako random ispadne tono 1.0

            // indeksi susjeda uz ciklicke rubne uvjete
            int iRight = (iX < lx - 1) ? iX + 1 : 0;
            int iDown = (iY < ly - 1) ? iY + 1 : 0;
            int iLeft = (iX > 0) ? iX - 1 : lx - 1;
            int iUp = (iY > 0) ? iY - 1 : ly - 1;

            // brojim kolko susjeda ima nizi stupac
            int neighbours = 0;
            if (grHeigthGrid[iRight, iY] < grHeigthGrid[iX, iY]) neighbours++;
            if (grHeigthGrid[iX, iDown] < grHeigthGrid[iX, iY]) neighbours++;
            if (grHeigthGrid[iLeft, iY] < grHeigthGrid[iX, iY]) neighbours++;
            if (grHeigthGrid[iX, iUp] < grHeigthGrid[iX, iY]) neighbours++;

            // nasumicno biram sudbinu kockice
            float rnd = Ranlux.NextFloat();

            // kristalizacija
            if (rnd < grTransitionProbability[5])
            {
                grHeigthGrid[iX, iY]++;
            }
            // evaporacija
            else if (rnd < grTransitionProbability[5] + grTransitionProbability[neighbours])
            {
                grHeigthGrid[iX, iY]--;
            }

        }

        // monte carlo loop
        public void grMonteCarloRun()
        {
            for (int j = 0; j < lx * ly; j++)
            {
                grMonteCarloStep();
            }

            grCalculateEnergy();
        }
	
        void UpdateGrain()
        {
			grMonteCarloRun();
			grTimeStepCount++;
        }



        double adtemperature = 1;
		const double adFluxMonolayerFrequency = 1.8e-9f;

		private const bool FLAG_grainTypeOlivin = false;

		//energy in  kelvin
	    const double adEnergyDiffuse = FLAG_grainTypeOlivin ? 287 : 511;
	    const double adEnergyDesorb = FLAG_grainTypeOlivin ? 373 : 658;

		/**** Adsorbcija vodika ***/

		int[,] adGrid;
		List<adatom> adList = new List<adatom>();

        double adDepositionWaitingTime = 0;
		double[] adDiffusionTime = new double[5];
		double[] adDesorptionTime = new double[5];
		double[] adActionWaitingTime = new double[5];
		double[] adFractionForDesorption = new double[5];
		double adLateralBondStrength = 0;
		
		long adTimeStepCount = 0;
		double adTimeStep = 0;
		int adHydrogenAtomsIncoming = 0;
		int adHydrogenAtomsDeposited = 0;
		int adHydrogenAtomsDesorped = 0;
		int adHydrogenMoleculeCount = 0;

		double mcTimespanAdatom;
		double elapsedTimeAdatom;
		double adTotalElapsedTime;

		public bool FLAG_VariableTimeStep = false;

		public double AdTemperature
		{
			get { return adtemperature; }
			set
			{
				adtemperature = value > 0 ? value : 0;
				adCalculateTransitionProbabilities();
			}
		}
		public int AdNumberOfInstances
        {
            get
            {
                return adList.Count;
            }
        }
		public double AdTimeStep
        {
            get
            {
				return adTimeStep;
            }
        }
		public double AdTotalElapsedTime
		{
			get { return adTotalElapsedTime; }
		}
		public long AdElapsedTimeSteps
		{
			get { return adTimeStepCount; }
		}
		public int AdHydrogenAtomsIncoming
		{
			get { return adHydrogenAtomsIncoming; }
		}
		public int AdHydrogenAtomsDeposited
		{
			get { return adHydrogenAtomsDeposited; }
		}
		public int AdHydrogenAtomsDesorped
		{
			get { return adHydrogenAtomsDesorped; }
		}
		public int AdHydrogenMoleculeCount
		{
			get { return adHydrogenMoleculeCount; }
		}

		public double AdTimeOfDeposition
		{
			get { return adDepositionWaitingTime; }
		}
		public double[] AdTimeOfAction
		{
			get { return adActionWaitingTime; }
		}
		public double AdLateralBondStrength
		{
			get { return adLateralBondStrength; }
			set
			{
				adLateralBondStrength = value >= 0 ? value : 0;
				adCalculateTransitionProbabilities();
			}
		}
		public double[] AdTimeOfDiffusion
		{
			get { return adDiffusionTime; }
		}
		public double[] AdTimeOfDesorption
		{
			get { return adDesorptionTime; }
		}
		public double[] AdFractionForDesorption
		{
			get { return adFractionForDesorption; }
		}
		

		class surfaceAtom
		{
			public int x = 0;
			public int y = 0;
			public double waitTime = 0; // seconds
			public double elapsedTime = 0; // seconds

			public surfaceAtom() { }
			public surfaceAtom(int x, int y, double waitTime, double elapsedTime)
			{
				this.x = x;
				this.y = y;
				this.waitTime = waitTime;
				this.elapsedTime = elapsedTime;
			}
		}
		class adatom : surfaceAtom
		{

			private static double[] adActionWaitingTime = new double[5];

			public adatom(double[] adActionWaitingTime, int neighbours)
			{
				adatom.adActionWaitingTime = adActionWaitingTime;

				NextWaitTime(neighbours);
			}
			public adatom(int x, int y, double[] adActionWaitingTime, int neighbours)
			{
				adatom.adActionWaitingTime = adActionWaitingTime;

				this.x = x;
				this.y = y;
				NextWaitTime(neighbours);
			}

			public void NextWaitTime(int neighbours)
			{
				elapsedTime = 0;
				waitTime = -adActionWaitingTime[neighbours] * Math.Log(Ranlux.NextDouble());
			}

			
		}

        // inicijalizira mrezu
        public void adInitialize()
        {
			adList = new List<adatom>();

            adGrid = new int[lx, ly];

            for (int i = 0; i < lx; i++)
            {
                for (int j = 0; j < ly; j++)
                {
                    adGrid[i, j] = 0;
                }
            }
			adTotalElapsedTime = 0;
			adTimeStepCount = 0;
			mcTimespanAdatom = 0;
			efficiency = 0;
			adHydrogenAtomsDeposited = 0;
			adHydrogenAtomsDesorped = 0;
			adHydrogenAtomsIncoming = 0;
			adHydrogenMoleculeCount = 0;
			FLAG_adCalculateNextDraw = true;

			InitializeAnim();
			EfficiencyInitialize();
			adCalculateTransitionProbabilities();
        }

        // inicijalizira frekvencije adsorpcije i desorpcije
        public void adCalculateTransitionProbabilities()
        {
			double adDepositionFrequency = lx * ly * adFluxMonolayerFrequency;
			adDepositionWaitingTime = 1.0 / adDepositionFrequency;

			//Console.WriteLine("adDepositionWaitingTime: " + adDepositionWaitingTime);

			for (int i = 0; i < 5; i++)
			{
				double LateralEnergy = i * adLateralBondStrength * adEnergyDesorb;
				double adDiffusionFrequency = 1e12 * Math.Exp(-(adEnergyDiffuse + LateralEnergy) / AdTemperature);
				double adDesorbtionFrequency = 1e12 * Math.Exp(-(adEnergyDesorb + LateralEnergy) / AdTemperature);

				

				adDiffusionTime[i] = 1.0 / adDiffusionFrequency;
				adDesorptionTime[i] = 1.0 / adDesorbtionFrequency;
				adActionWaitingTime[i] = 1.0 / (adDiffusionFrequency + adDesorbtionFrequency);
				adFractionForDesorption[i] = adDesorbtionFrequency * adActionWaitingTime[i];

				/*
				
				Console.WriteLine("adDiffusionTime[" + i + "]: " + adDiffusionTime[i]);
				Console.WriteLine("adDesorptionTime[" + i + "]: " + adDesorptionTime[i]);
				Console.WriteLine("adActionWaitingTime[" + i + "]: " + adActionWaitingTime[i]);
				Console.WriteLine("adLimitForDesorption[" + i + "]: " + adFractionForDesorption[i]);
				*/
			}
        }
		
		public int adHighNeighbourNumber(int iX, int iY)
		{
			int iRight = (iX < lx - 1) ? iX + 1 : 0;
			int iDown = (iY < ly - 1) ? iY + 1 : 0;
			int iLeft = (iX > 0) ? iX - 1 : lx - 1;
			int iUp = (iY > 0) ? iY - 1 : ly - 1;

			// brojim kolko ima visih susjeda
			int neighbours = 0;
			if (grHeigthGrid[iRight, iY] > grHeigthGrid[iX, iY]) neighbours++;
			if (grHeigthGrid[iX, iDown] > grHeigthGrid[iX, iY]) neighbours++;
			if (grHeigthGrid[iLeft, iY] > grHeigthGrid[iX, iY]) neighbours++;
			if (grHeigthGrid[iX, iUp] > grHeigthGrid[iX, iY]) neighbours++;

			return neighbours;
		}

		// pomakni vrijeme i updejtaj adatome
		void UpdateAdatom()
        {
			updateTimeStep();

			double timeStep = AdTimeStep;

			// inkrement vremena
			elapsedTimeAdatom += timeStep;
			adTotalElapsedTime += timeStep;
			
            if (elapsedTimeAdatom > mcTimespanAdatom)
            {
                int x = (int)(Ranlux.NextFloat() * lx);
                int y = (int)(Ranlux.NextFloat() * ly);

                if (adGrid[x, y] == 0)
                {
					if (FLAG_desorbRun)
					{
						HadsList.Add(new surfaceAtom(x, y, 10, 0));
					}
					else
					{
						adList.Add(new adatom(x, y, adActionWaitingTime, adHighNeighbourNumber(x, y)));
					}

                    adGrid[x, y]++;
					adHydrogenAtomsDeposited++;
                }
				adHydrogenAtomsIncoming++;
                elapsedTimeAdatom = 0;
                mcTimespanAdatom = -adDepositionWaitingTime * Math.Log(Ranlux.NextDouble());
            }


            // update adatom time
            for (int i = 0; i < adList.Count; i++)
            {
				adList[i].elapsedTime += timeStep;
            }

			// ako treba skenirati cijelu mrezu za naci dva atoma na istom mjestu
			bool FLAG_removeMolecule = false;

            // perform adatom action
			for (int i = adList.Count - 1; i >= 0 && adList.Count > 0; i--)
            {
                if (adList[i].elapsedTime > adList[i].waitTime)
                {
                    double decideDesorp = Ranlux.NextDouble();

                    // desorb
					if (decideDesorp < adFractionForDesorption[adHighNeighbourNumber(adList[i].x, adList[i].y)])
                    {
						adGrid[adList[i].x, adList[i].y]=0;

						if (FLAG_desorbRun)
						{
							HdesList.Add(new surfaceAtom(adList[i].x, adList[i].y, 30, 0));
						}

                        adList.RemoveAt(i);
						adHydrogenAtomsDesorped++;

                    }
                    else //move
                    {
                        float decideMove = Ranlux.NextFloat();

                        // makni atom sa starog mjesta
                        adGrid[adList[i].x, adList[i].y]--;

                        if (decideMove < 0.25f)
                        {
                            adList[i].x = (adList[i].x < lx - 1) ? adList[i].x + 1 : 0;
                        }
                        else if (decideMove < 0.5f)
                        {
                            adList[i].x = (adList[i].x > 0) ? adList[i].x - 1 : lx - 1;
                        }
                        else if (decideMove < 0.75f)
                        {
                            adList[i].y = (adList[i].y < ly - 1) ? adList[i].y + 1 : 0;
                        }
                        else
                        {
                            adList[i].y = (adList[i].y > 0) ? adList[i].y - 1 : ly - 1;
                        }

                        adGrid[adList[i].x, adList[i].y]++;

						adList[i].NextWaitTime(adHighNeighbourNumber(adList[i].x, adList[i].y));

                        // ako se nadju dva na istom mjestu
                        if (adGrid[adList[i].x, adList[i].y] > 1)
                        {
                            FLAG_removeMolecule = true;
                        }
                    }

                }
            }
            // remove molecules
            if (FLAG_removeMolecule)
            {
                for (int i = 0; i < lx; i++)
                {
                    for (int j = 0; j < ly; j++)
                    {
                        if (adGrid[i, j] >= 2)
                        {
							adList.RemoveAll(adatom => adatom.x == i && adatom.y == j);
                            adGrid[i, j] = 0;
                            adHydrogenMoleculeCount++;

							if (FLAG_desorbRun)
							{
								H2List.Add(new surfaceAtom(i, j, 30, 0));
							}
                        }
                    }
                }
            }

			adTimeStepCount++;

        }
		void updateTimeStep()
		{
			if (FLAG_VariableTimeStep)
			{
				double min = double.MaxValue;

				for (int i = 0; i < adList.Count; i++)
				{
					if (adList[i].waitTime < min) min = adList[i].waitTime;
					//if (adList[i].waitTime < min) min = adList[i].waitTime - adList[i].elapsedTime; // mozda bolje ali pazi da nije 0
				}


				adTimeStep = Math.Min(adDepositionWaitingTime, min);
			}
			else
			{
				adTimeStep = Math.Min(adDepositionWaitingTime, Math.Min(adDiffusionTime[0], adDesorptionTime[0])) / 2;
			}
		}


		/**** Animacije ****/

		List<surfaceAtom> H2List = new List<surfaceAtom>();
		List<surfaceAtom> HdesList = new List<surfaceAtom>();
		List<surfaceAtom> HadsList = new List<surfaceAtom>();
		void InitializeAnim()
		{
			H2List = new List<surfaceAtom>();
			HdesList = new List<surfaceAtom>();
			HadsList = new List<surfaceAtom>();
		}
		
		void UpdateAnim()
		{
			for (int i = H2List.Count-1; i >= 0; i--)
			{
				H2List[i].elapsedTime++;

				if (H2List[i].elapsedTime >= H2List[i].waitTime)
				{
					H2List.RemoveAt(i);
				}
			}
			for (int i = HdesList.Count - 1; i >= 0; i--)
			{
				HdesList[i].elapsedTime++;

				if (HdesList[i].elapsedTime >= HdesList[i].waitTime)
				{
					HdesList.RemoveAt(i);
				}
			}
			for (int i = HadsList.Count - 1; i >= 0; i--)
			{
				HadsList[i].elapsedTime++;

				if (HadsList[i].elapsedTime >= HadsList[i].waitTime)
				{
					adList.Add(new adatom(HadsList[i].x, HadsList[i].y, adActionWaitingTime, adHighNeighbourNumber(HadsList[i].x, HadsList[i].y)));
					HadsList.RemoveAt(i);
				}
			}

		}

		/*** Ucinkovitost ***/

		//int efDeposited = 0;
		int efCycleFormedMoleculesPrev = 0;
		int efCycleIncomingHydrogenPrev = 0;
		int efCycleDepositedPrev = 0;
		int efCycleTrigger = 1000;
		int efCycle = 0;

		int efCycleAdatomsPrev = 0;
		double efficiency = 0;
		public bool FLAG_steadyState = false;
		public bool FLAG_meanEfficiencyFinish = false;
	    int adEfficiencyPassCounter
	    {
			get
			{
				if (LengthX*LengthY < 1000)
				{
					return 100;
				}
				else if (LengthX*LengthY < 10000)
				{
					return 30;
				}
				else
				{
					return 15;
				}
				
			}
	    }
		int efAveragingPass = 0;
		double meanEfficiency = 0;
		

		public double AdEfficiency
		{
			get { return efficiency; }
		}
		public double AdMeanEfficiency
		{
			get { return meanEfficiency; }
		}
		public int AdEfCycleInDepositions
		{
			get
			{
				//return efCycleTrigger;
				return (int) Math.Max((LengthX*LengthY/100f), 10);
			}
		}
		void EfficiencyInitialize()
		{
			FLAG_meanEfficiencyFinish = false;
			FLAG_steadyState = false;
			efCycleDepositedPrev = 0;
			efCycleAdatomsPrev = AdInstancesCount;
			efCycleIncomingHydrogenPrev = AdHydrogenAtomsIncoming;
			efCycleFormedMoleculesPrev = AdHydrogenMoleculeCount;
			efficiency = 0;
			efAveragingPass = 0;
			efCycle = 0;
			meanEfficiency = 0;
		}

		void UpdateEfficiency()
		{
			int counter = AdHydrogenAtomsIncoming - efCycleIncomingHydrogenPrev;

			if (counter == AdEfCycleInDepositions)
			{
				efCycle++;

				int cycleIncoming = AdHydrogenAtomsIncoming - efCycleIncomingHydrogenPrev;
				int cycleFormed = AdHydrogenMoleculeCount - efCycleFormedMoleculesPrev;

				if (cycleIncoming != 0)
				{
					efficiency = 2 * (double)cycleFormed / cycleIncoming;
				}
				else
				{
					efficiency = 0;
				}

				if (!FLAG_steadyState && Math.Abs((AdInstancesCount - (float)efCycleAdatomsPrev) / (float)efCycleAdatomsPrev) < 0.01f)
				{
					FLAG_steadyState = true;
				}

				/*
				// ako je efikasnost nakon 10 ciklusa jos ucijek nula
				if (efCycle > 10 && efficiency < 0.001)
				{
					meanEfficiency = 0;
					FLAG_meanEfficiencyFinish = true;
				}
				*/

				// biljezi za konacni rezultat
				if (FLAG_steadyState && !FLAG_meanEfficiencyFinish)
				{
					efAveragingPass++;
					meanEfficiency += efficiency;
					//Console.WriteLine(meanEfficiency);

					if (efAveragingPass == adEfficiencyPassCounter)
					{
						meanEfficiency /= adEfficiencyPassCounter;
						efAveragingPass = 0;
						FLAG_meanEfficiencyFinish = true;
					}
				}

				efCycleAdatomsPrev = AdInstancesCount;
				efCycleDepositedPrev = AdHydrogenAtomsDeposited;
				efCycleIncomingHydrogenPrev = AdHydrogenAtomsIncoming;
				efCycleFormedMoleculesPrev = AdHydrogenMoleculeCount;
				counter = 0;
			}

			
		}

		double lowerTemp=FLAG_grainTypeOlivin ? 3 : 8;
		double upperTemp=FLAG_grainTypeOlivin ? 80 : 100;


		double increment=1;
		public bool FLAG_TemperatureRangeRun = false;

		public void InitializeTemperatureRange()
		{
			StreamWriter sw = new StreamWriter("eff.txt", false);
			sw.WriteLine("Temperature\tEfficiency\tAdatoms");
			sw.Close();


			FLAG_TemperatureRangeRun = true;

			AdTemperature = lowerTemp;
		}
		public void TemperatureRangeUpdate()
		{
			//FLAG_RangeRun = true;
			StreamWriter sw = new StreamWriter("eff.txt", true);
			sw.WriteLine(AdTemperature + "\t" + AdMeanEfficiency + "\t" + AdInstancesCount);
			sw.Close();

			if (AdTemperature + increment > upperTemp)
			{
				FLAG_TemperatureRangeRun = false;
				bw.CancelAsync();
				//Application.Exit();
				return;
			}

			AdTemperature += increment;

			//sw.Close();
			
			

			adInitialize();
		}

		public bool FLAG_SizeRangeRun = false;
		public void InitializeSizeRange()
		{
			StreamWriter sw = new StreamWriter("eff-size.txt", true);
			sw.WriteLine("N\tEfficiency\t");
			sw.Close();

			FLAG_SizeRangeRun = true;

			LengthX = 5;
			LengthY = 5;
		}
		public void SizeRangeUpdate()
		{
			//FLAG_RangeRun = true;
			StreamWriter sw = new StreamWriter("eff-size.txt", true);
			sw.WriteLine(LengthX*LengthY + "\t" + AdMeanEfficiency);
			sw.Close();

			if (LengthX * LengthY > 500*500)
			{
				FLAG_SizeRangeRun = false;
				bw.CancelAsync();
				//Application.Exit();
				return;
			}

			if (LengthX < 100)
			{
				LengthX += 1;
				LengthY += 1;
			}
			else
			{
				LengthX += 10;
				LengthY += 10;
			}
			//AdTemperature += increment;
			
	
			
			//LengthX = (int) (LengthX*1.5f);
			//LengthY = (int) (LengthY*1.5f);

	
			adInitialize();
		}

		public bool FLAG_TimeRangeRun = false;
	    private int timerecordingStep = 0;
		public void InitializeTimeRange()
		{
			StreamWriter sw = new StreamWriter("time.txt", false);
			sw.WriteLine("Time\tAdatoms\tMoleculesFormed\tSteadyStep");
			sw.Close();

			timerecordingStep = int.MaxValue-1;
			FLAG_TimeRangeRun = true;
			
			adInitialize();
		}
		public void TimeRangeUpdate()
		{
			timerecordingStep++;
			if (timerecordingStep >= 100)
			{
				StreamWriter sw = new StreamWriter("time.txt", true);
				sw.WriteLine(AdTotalElapsedTime + "\t" + AdInstancesCount + "\t" + AdHydrogenMoleculeCount + "\t" + (FLAG_steadyState ? 1:0));
				sw.Close();

				if (FLAG_meanEfficiencyFinish)
				{
					//FLAG_TimeRangeRun = false;
					//bw.CancelAsync();
				}

				timerecordingStep = 0;
			}
			
		}

		double latStartTemp = FLAG_grainTypeOlivin ? 8 : 14;
		public bool FLAG_LateralRangeRun = false;
	    private double lastEfficiency = 0;
		public void InitializeLateralRange()
		{

			StreamWriter sw = new StreamWriter("lateral.txt", false);

			sw.WriteLine("Lateral\tTemperature\tEfficiency");
			sw.Close();

			FLAG_LateralRangeRun = true;
			lastEfficiency = 0;
			adLateralBondStrength = 0;
			AdTemperature = latStartTemp;
		}
		public void LateralRangeUpdate()
		{
			if (lastEfficiency >= 0.5 && AdMeanEfficiency < 0.5)
			{
				StreamWriter sw = new StreamWriter("lateral.txt", true);

				double e1, e2, t1, t2, t50;

				e1 = lastEfficiency;
				e2 = AdMeanEfficiency;
				t1 = AdTemperature - 1;
				t2 = AdTemperature;

				t50 = t2 - (0.5 - e2)/(e1 - e2)*(t2 - t1);

				sw.WriteLine(adLateralBondStrength + "\t" + t50 + "\t" + lastEfficiency);
				sw.Close();

				if (adLateralBondStrength >= 0.99)
				{
					FLAG_LateralRangeRun = false;
					bw.CancelAsync();
					//Application.Exit();
					return;
				}
				adLateralBondStrength += 0.1;
				AdTemperature = latStartTemp;
				lastEfficiency = 0;
			}
			else
			{
				lastEfficiency = AdMeanEfficiency;
				AdTemperature += 1;
			}

			adInitialize();

		}


		/*** Crtanje ***/

		public bool FLAG_grainRun = false;
		public bool FLAG_adatomRun = false;
		public bool FLAG_desorbRun = true;
		public bool FLAG_grCalculateNextDraw = true;
		public bool FLAG_adCalculateNextDraw = true;
		//private bool FLAG_h2CalculateNextDraw = true;
		public bool FLAG_Draw = true;

		public int GrInstancesCount
		{
			get
			{
				return instanceTransformsGrain != null ? instanceTransformsGrain.Length : 0;
			}
		}
		public int AdInstancesCount
		{
			get
			{
				return adList.Count;
			}
		}

		void UpdateAll()
		{
			if (FLAG_grainRun)
			{
				UpdateGrain();
				FLAG_grCalculateNextDraw = true;
			}

			if (FLAG_adatomRun)
			{
				UpdateAdatom();
				FLAG_adCalculateNextDraw = true;
			}

			if (FLAG_desorbRun)
			{
				UpdateAnim();
				//FLAG_h2CalculateNextDraw = true;
			}
			else
			{
				if (H2List.Count > 0) H2List.Clear();
				if (HdesList.Count > 0) HdesList.Clear();
				if (HadsList.Count > 0) HadsList.Clear();
			}

			UpdateEfficiency();

			if (FLAG_TemperatureRangeRun)
			{
				if (FLAG_meanEfficiencyFinish)
				{
					TemperatureRangeUpdate();
				}
			}
			if (FLAG_SizeRangeRun)
			{
				if (FLAG_meanEfficiencyFinish)
				{
					SizeRangeUpdate();
				}
			}
			if (FLAG_TimeRangeRun)
			{
				TimeRangeUpdate();	
			}
			if (FLAG_LateralRangeRun)
			{
				if (FLAG_meanEfficiencyFinish)
				{
					LateralRangeUpdate();
				}

			}
		}

		// popuni matrice pozicija atoma za iscrtavanje
		void calculateGrainInstanceTransforms()
		{
			int padding = 0;

			int count = 0;
			//couting pass
			for (int i = 0; i < lx; i++)
			{
				for (int j = 0; j < ly; j++)
				{
					count++;
					for (int k = -padding; k < grHeigthGrid[i, j]; k++)
					{
						count++;
					}
				}
			}
			//// Gather instance transform matrices into a single array.
			//Array.Resize(ref instanceTransforms, count);

			instanceTransformsGrain = new Matrix[count];

			count = 0;
			//configUpdate();
			for (int i = 0; i < lx; i++)
			{
				for (int j = 0; j < ly; j++)
				{
					//DrawModel(kocka, Matrix.CreateTranslation(new Vector3(i, ih[i,j], j)));
					instanceTransformsGrain[count++] = Matrix.CreateTranslation(new Vector3(i, grHeigthGrid[i, j], j));
					for (int k = -padding; k < grHeigthGrid[i, j]; k++)
					{
						//DrawModel(kocka, Matrix.CreateTranslation(new Vector3(i, k, j)));
						instanceTransformsGrain[count++] = Matrix.CreateTranslation(new Vector3(i, k, j));
					}
				}
			}
		}
		void calculateAdatomInstanceTransforms()
		{

			instanceTransformsAd = new Matrix[adList.Count];

			for (int k = 0; k < adList.Count; k++)
			{
				instanceTransformsAd[k] = Matrix.CreateTranslation(new Vector3(
					adList[k].x, grHeigthGrid[adList[k].x, adList[k].y] + 1, adList[k].y));

			}
		}

		const float finalDesorbHeigth = 30;
		const float finalAdsorbHeigth = 20;
		void calculateH2InstanceTransforms()
		{
			instanceTransformsH2 = new Matrix[H2List.Count];

			for (int k = 0; k < H2List.Count; k++)
			{
				float completed = (float)(H2List[k].elapsedTime / H2List[k].waitTime);

				instanceTransformsH2[k] = Matrix.CreateScale(1 - 0.6f*completed * completed) * Matrix.CreateTranslation(new Vector3(
					H2List[k].x, grHeigthGrid[H2List[k].x, H2List[k].y] + 1 + finalDesorbHeigth * completed, H2List[k].y));
			}
		}
		void calculateHdesInstanceTransforms()
		{
			instanceTransformsHdes = new Matrix[HdesList.Count];

			for (int k = 0; k < HdesList.Count; k++)
			{
				float completed = (float)(HdesList[k].elapsedTime / HdesList[k].waitTime);

				instanceTransformsHdes[k] = Matrix.CreateScale(1 -  0.6f*completed * completed) * Matrix.CreateTranslation(new Vector3(
					HdesList[k].x, grHeigthGrid[HdesList[k].x, HdesList[k].y] + 1 + finalDesorbHeigth * completed, HdesList[k].y));
			}
		}
		void calculateHadsInstanceTransforms()
		{
			instanceTransformsHads = new Matrix[HadsList.Count];

			for (int k = 0; k < HadsList.Count; k++)
			{
				float completed = (float)(HadsList[k].elapsedTime / HadsList[k].waitTime);

				instanceTransformsHads[k] = Matrix.CreateScale(0.4f + 0.6f*completed * completed) * Matrix.CreateTranslation(new Vector3(
					HadsList[k].x, grHeigthGrid[HadsList[k].x, HadsList[k].y] + 1 + finalAdsorbHeigth * (1 - completed), HadsList[k].y));
			}
		}

        protected override void Draw()
        {
			DepthStencilState depthStencilState = new DepthStencilState();
			depthStencilState.DepthBufferEnable = true;
			GraphicsDevice.DepthStencilState = depthStencilState;


			//GraphicsDevice.

            GraphicsDevice.Clear(Color.Black );
			//GraphicsDevice.Clear(Color.White);

            GraphicsDevice.DepthStencilState = DepthStencilState.Default;
            GraphicsDevice.BlendState = BlendState.Opaque;			

			if (FLAG_Draw)
			{
				if (FLAG_grCalculateNextDraw)
				{
					calculateGrainInstanceTransforms();
					FLAG_grCalculateNextDraw = false;
				}
				if (FLAG_adCalculateNextDraw)
				{
					calculateAdatomInstanceTransforms();
					FLAG_adCalculateNextDraw = false;
				}
				calculateH2InstanceTransforms();
				calculateHdesInstanceTransforms();
				calculateHadsInstanceTransforms();

				DrawModelHardwareInstancing(instancedModel_grain, cubeColor.Grain, instancedModelBones_grain, instanceTransformsGrain, camera.view, camera.projection);
				DrawModelHardwareInstancing(instancedModel_ad, cubeColor.Hydrogen, instancedModelBones_ad, instanceTransformsHads, camera.view, camera.projection);
				DrawModelHardwareInstancing(instancedModel_ad, cubeColor.Hydrogen, instancedModelBones_ad, instanceTransformsAd, camera.view, camera.projection);
				DrawModelHardwareInstancing(instancedModel_h2, cubeColor.Red, instancedModelBones_h2, instanceTransformsH2, camera.view, camera.projection);
				DrawModelHardwareInstancing(instancedModel_ad, cubeColor.Hydrogen, instancedModelBones_ad, instanceTransformsHdes, camera.view, camera.projection);


			}

        }

		protected override void OnMouseEnter(System.EventArgs e)
		{
			this.Focus();
			base.OnMouseEnter(e);
		}
		protected override void OnMouseMove(MouseEventArgs e)
		{
			camera.Update(e);
			base.OnMouseMove(e);
		}
		protected override void OnMouseWheel(MouseEventArgs e)
		{
			camera.Update(e);
			base.OnMouseWheel(e);
		}

		private Model LoadModel(string assetName, out Texture2D[] textures)
		{
			Model newModel = Content.Load<Model>(assetName);
			
			int textureNum = 0;
			foreach (ModelMesh mesh in newModel.Meshes)
				foreach (BasicEffect currentEffect in mesh.Effects)
					textureNum++;

			textures = new Texture2D[textureNum];
			int i = 0;
			foreach (ModelMesh mesh in newModel.Meshes)
				foreach (BasicEffect currentEffect in mesh.Effects)
					textures[i++] = currentEffect.Texture;

			/*
			foreach (ModelMesh mesh in newModel.Meshes)
				foreach (ModelMeshPart meshPart in mesh.MeshParts)
					meshPart.Effect = basicEffect.Clone();
			*/
			return newModel;
		}
		private void DrawModel(Model model, Matrix wMatrix)
		{
			SamplerState ss = new SamplerState();
			ss.Filter = TextureFilter.Anisotropic;
			ss.MaxAnisotropy = 16;
			ss.AddressU = TextureAddressMode.Clamp;
			ss.AddressV = TextureAddressMode.Clamp;

			GraphicsDevice.SamplerStates[0] = ss;

			Matrix[] modelTransforms = new Matrix[model.Bones.Count];
			model.CopyAbsoluteBoneTransformsTo(modelTransforms);
			//int i = 0;
			foreach (ModelMesh mesh in model.Meshes)
			{
				foreach (BasicEffect basicEffect in mesh.Effects)
				{
					basicEffect.View = camera.view;
					basicEffect.Projection = camera.projection;
					basicEffect.World = modelTransforms[mesh.ParentBone.Index] * wMatrix;

					//basicEffect.VertexColorEnabled = true;
					basicEffect.EnableDefaultLighting();
					//Console.WriteLine(modelTextures.Length);
					//basicEffect.Texture = textures[i++];
					//basicEffect.TextureEnabled = true;

				}
				mesh.Draw();
			}
		}

        Model instancedModel_grain;
        Matrix[] instancedModelBones_grain;
		Model instancedModel_ad;
		Matrix[] instancedModelBones_ad;
		Model instancedModel_h2;
		Matrix[] instancedModelBones_h2;

        Effect instanceModelEffect;
        Matrix[] instanceTransformsGrain;
		Matrix[] instanceTransformsAd;
		Matrix[] instanceTransformsH2;
		Matrix[] instanceTransformsHdes;
		Matrix[] instanceTransformsHads;
        DynamicVertexBuffer instanceVertexBuffer;

        // To store instance transform matrices in a vertex buffer, we use this custom
        // vertex type which encodes 4x4 matrices as a set of four Vector4 values.
        static VertexDeclaration instanceTransformVertexDeclaration = new VertexDeclaration(
            new VertexElement(0, VertexElementFormat.Vector4, VertexElementUsage.BlendWeight, 0),
            new VertexElement(16, VertexElementFormat.Vector4, VertexElementUsage.BlendWeight, 1),
            new VertexElement(32, VertexElementFormat.Vector4, VertexElementUsage.BlendWeight, 2),
            new VertexElement(48, VertexElementFormat.Vector4, VertexElementUsage.BlendWeight, 3)
        );

        static class cubeColor
        {
            public static Vector4 Red = new Vector4(1, 0, 0, 1);
            public static Vector4 Green = new Vector4(0, 1, 0, 1);
            public static Vector4 Blue = new Vector4(0, 0, 1, 1);

            public static Vector4 Grain = new Vector4(0, 0.4f, 0.3f, 1);
            public static Vector4 Hydrogen = new Vector4(0.8f, 1f, 1f, 1);

			//sivo
			//public static Vector4 Grain = new Vector4(1f, 1f, 1f, 1);
			//public static Vector4 Hydrogen = new Vector4(0.3f, 0.3f, 1f, 1);

			public static Vector4 FromColor(Color color)
			{
				return new Vector4((float)color.R / byte.MaxValue,
					(float)color.G / byte.MaxValue, 
					(float)color.B / byte.MaxValue,
					(float)color.A / byte.MaxValue);
			}

        }
		/// <summary>
		/// Efficiently draws several copies of a piece of geometry using hardware instancing.
		/// </summary>
        void DrawModelHardwareInstancing(Model model, Vector4 color, Matrix[] modelBones, Matrix[] instances, Matrix view, Matrix projection)
        {
			if (instances== null || instances.Length == 0)
                return;

            // If we have more instances than room in our vertex buffer, grow it to the neccessary size.
            if ((instanceVertexBuffer == null) ||
                (instances.Length > instanceVertexBuffer.VertexCount))
            {
                if (instanceVertexBuffer != null)
                    instanceVertexBuffer.Dispose();


                instanceVertexBuffer = new DynamicVertexBuffer(GraphicsDevice, instanceTransformVertexDeclaration,
                                                               instances.Length, BufferUsage.WriteOnly);
            }

            // Transfer the latest instance transform matrices into the instanceVertexBuffer.
            instanceVertexBuffer.SetData(instances, 0, instances.Length, SetDataOptions.Discard);

            foreach (ModelMesh mesh in model.Meshes)
            {
                foreach (ModelMeshPart meshPart in mesh.MeshParts)
                {
                    // Tell the GPU to read from both the model vertex buffer plus our instanceVertexBuffer.
                    GraphicsDevice.SetVertexBuffers(
                        new VertexBufferBinding(meshPart.VertexBuffer, meshPart.VertexOffset, 0),
                        new VertexBufferBinding(instanceVertexBuffer, 0, 1)
                    );

                    GraphicsDevice.Indices = meshPart.IndexBuffer;
					
                    // Set up the instance rendering effect.
                    //Effect effect = meshPart.Effect;

                    Effect effect = instanceModelEffect;

					//RasterizerState rs = new RasterizerState();
					//rs.CullMode = CullMode.CullCounterClockwiseFace;
					//GraphicsDevice.RasterizerState = rs;

                    effect.Parameters["VertexColor"].SetValue(color);


					if (color == cubeColor.Grain)
					{
						effect.CurrentTechnique = effect.Techniques["HardwareInstancing2"];
					}
					else
					{
						effect.CurrentTechnique = effect.Techniques["HardwareInstancing"];
					}
                    
                    effect.Parameters["World"].SetValue(modelBones[mesh.ParentBone.Index]);
                    effect.Parameters["View"].SetValue(view);
                    effect.Parameters["Projection"].SetValue(projection);

                    // Draw all the instance copies in a single call.
                    foreach (EffectPass pass in effect.CurrentTechnique.Passes)
                    {
                        pass.Apply();

                        GraphicsDevice.DrawInstancedPrimitives(PrimitiveType.TriangleList, 0, 0,
                                                               meshPart.NumVertices, meshPart.StartIndex,
                                                               meshPart.PrimitiveCount, instances.Length);
                    }
                }
            }
        }

    }
}
