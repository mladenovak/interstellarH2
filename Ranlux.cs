using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;

namespace InterstellarIce
{
	static class Ranlux
	{
		const string dllPath = "Ranlux.dll";

		[DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
		static extern void ranlxs(float[] r, int n);
		[DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
		static extern void rlxs_init(int level, int seed);
		[DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
		static extern int rlxs_size();
		[DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
		static extern void rlxs_get(int[] state);
		[DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
		static extern void rlxs_reset(int[] state);

		[DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
		static extern void ranlxd(double[] r, int n);
		[DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
		static extern void rlxd_init(int level, int seed);
		[DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
		static extern int rlxd_size();
		[DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
		static extern void rlxd_get(int[] state);
		[DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
		static extern void rlxd_reset(int[] state);

		[DllImport(dllPath, CallingConvention = CallingConvention.Cdecl)]
		static extern int main();

		/// <summary>
		/// 
		/// </summary>
		/// <param name="level">Luxury level: 0, 1 or 2</param>
		/// <param name="seed">Integer between 1 and 2^31</param>
		public static void InitializeFloatGenerator(int level, int seed)
		{
			rlxs_init(level, seed);
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="level">Luxury level: 1 or 2</param>
		/// <param name="seed">Integer between 1 and 2^31</param>
		public static void InitializeDoubleGenerator(int level, int seed)
		{
			rlxd_init(level, seed);
		}

        // otvoreni itnerval: nema 0 ili 1
		public static double NextDouble()
		{
			double[] r = new double[1];

            do
            {
                ranlxd(r, 1);
            } while (r[0] == 0);

			return r[0];
		}
        public static float NextFloat()
        {
            float[] r = new float[1];

            do
            {
                ranlxs(r, 1);
            } while (r[0] == 0);

            return r[0];
        }

	};
}
