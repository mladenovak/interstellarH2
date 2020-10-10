using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using Microsoft.Xna.Framework;
//using Microsoft.Xna.Framework.Input;
using Microsoft.Xna.Framework.Graphics;

using System.Windows.Forms;

namespace InterstellarIce
{

	public class WinFormCamera
	{
		public Matrix view { get; protected set; }
		public Matrix projection { get; protected set; }
		public float camera_fi { get; set; }
		public float camera_theta { get; set; }
		public float camera_gamma { get; set; }
		public float camera_radius { get; set; }
		public Vector3 camera_target { get; set; }

		public bool FLAG_cameraNeedRecalculating = true;

		Vector3 position = new Vector3(0, 0, 1);
		Vector3 direction = new Vector3(0, 0, -1);
		Vector3 up = new Vector3(0, 1, 0);

		//float camera_nearPlane = 0.1f;
		//float camera_farPlane = 1000f;
		float camera_nearPlane = 2f;
		float camera_farPlane = 2048f;

		GraphicsDevice graphicsDevice;

		MouseEventArgs ePrev;

		public WinFormCamera(GraphicsDevice graphicsDevice,
			float camera_fi,
			float camera_theta,
			float camera_gamma,
			float camera_radius,
			Vector3 camera_target)
		{
			this.camera_fi = camera_fi;
			this.camera_theta = camera_theta;
			this.camera_gamma = camera_gamma;
			this.camera_radius = camera_radius;
			this.camera_target = camera_target;
			this.graphicsDevice = graphicsDevice;

			Update(new MouseEventArgs(MouseButtons.None,0,0,0,0));
		}

		public void UpdatePerspective()
		{
			projection = Matrix.CreatePerspectiveFieldOfView(
				MathHelper.ToRadians(30),
				graphicsDevice.Viewport.AspectRatio,
				camera_nearPlane,
				camera_farPlane);
		}

		public void Update(MouseEventArgs e)
		{
			UpdatePerspective();

			// rotate camera
			if (e.Button == MouseButtons.Middle || e.Button == MouseButtons.Right)
			{
				camera_fi += (ePrev.X - e.X) / 100f;
				camera_theta += (ePrev.Y - e.Y) / 100f;

				FLAG_cameraNeedRecalculating = true;
			}
			
			// zoom camera
			int scroll = e.Delta;
			
			
			if (scroll != 0)
			{
				camera_radius -= scroll * camera_radius / 1000f;

				// zumiraj prema kursoru
				Vector3 cameraTargetOnScreen = graphicsDevice.Viewport.Project(
					camera_target, projection, view, Matrix.Identity);

				Vector3 mouseTarget = graphicsDevice.Viewport.Unproject(
					new Vector3(e.X, e.Y, cameraTargetOnScreen.Z), projection, view, Matrix.Identity);

				Vector3 mouseCameraTarget = mouseTarget - camera_target;

				camera_target += Math.Sign(scroll) * mouseCameraTarget / 10f;

				FLAG_cameraNeedRecalculating = true;
			}

			if (FLAG_cameraNeedRecalculating)
			{
				position = new Vector3(0, 0, 1);
				direction = new Vector3(0, 0, -1);
				up = new Vector3(0, 1, 0);

				Matrix rot = Matrix.CreateFromYawPitchRoll(camera_fi, camera_theta, camera_gamma);

				view = Matrix.CreateLookAt(
					Vector3.Transform(position, rot) * camera_radius + camera_target,
					camera_target,
					Vector3.Transform(up, rot));

				FLAG_cameraNeedRecalculating = false;
			}
			ePrev = e;
		}
	}
}
