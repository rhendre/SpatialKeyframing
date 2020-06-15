
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;
using System;
using System.Globalization;
using System.IO;
using UnityEngine;



public class InterpolationScript : MonoBehaviour
    {
    // Create Matrix from spatial keyframing.txt and current cursor position.
       System.Numerics.Quaternion rotation;

    UnityEngine.Vector3 cursorPosition = new UnityEngine.Vector3();
        public UnityEngine.Vector3 scale = new UnityEngine.Vector3(1, 1, 1);
        double[,] interpolatedMatrix = new double[3, 3];
        int r = 0;

    UnityEngine.Vector3[] cursors = new UnityEngine.Vector3[3];
    System.Numerics.Quaternion[] leftLegRotations = new System.Numerics.Quaternion[3];
    System.Numerics.Quaternion[] rightLegRotations = new System.Numerics.Quaternion[3];
    System.Numerics.Quaternion[] neckRotations = new System.Numerics.Quaternion[3];
    System.Numerics.Quaternion[] leftShoulderRotations = new System.Numerics.Quaternion[3];
    System.Numerics.Quaternion[] rightShoulderRotations = new System.Numerics.Quaternion[3];
    int count = 0;

    public void Start()
    {
        //read the spatial keyframes from the text file.
        string path =
             Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);
        string file_path = Path.Combine(path, "SpatialKeyFrame.txt");
        StreamReader readerStream = new StreamReader(file_path);
        transform.hasChanged = false;
        while (!readerStream.EndOfStream)
        {
            string line = readerStream.ReadLine();

            if (line != null)
            {


                //load the cursor positions in the array cursor.
                if (line.Contains("Cursor"))
                {
                    string[] str = line.Split(',');
                    float positionx, positiony, positionz = 0.0f;
                    float.TryParse(str[0].Replace("Cursor[", ""), System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positionx);
                    float.TryParse(str[1], System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positiony);
                    float.TryParse(str[2].Replace("]", ""), System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positionz);
                    Vector3 position = new Vector3(positionx, positiony, positionz);
                    GameObject keyFrame = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                    keyFrame.transform.localScale = new Vector3(0.1f, 0.1f, 0.1f);
                    keyFrame.gameObject.name = "Cursor" + count;

                    keyFrame.gameObject.transform.position = position;
                    cursors.SetValue(position, count);

                    count++;
                }

                else if( line.Contains("LeftLeg"))
                {

                    float positionx, positiony, positionz, scalarw = 0.0f;
                    string[] str = line.Split(',');
                    float.TryParse(str[0].Replace("LeftLeg[", ""), System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positionx);
                    float.TryParse(str[1], System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positiony);
                    float.TryParse(str[2], System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positionz);
                    float.TryParse(str[3].Replace("]", ""), System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out scalarw);

                    System.Numerics.Quaternion quat = new System.Numerics.Quaternion();
                    quat.X = positionx;
                    quat.Y = positiony;
                    quat.Z = positionz;
                    quat.W = scalarw;


                    leftLegRotations[count] = quat;
                 
                }
                else if (line.Contains("RightLeg"))
                {

                    float positionx, positiony, positionz, scalarw = 0.0f;
                    string[] str = line.Split(',');
                    float.TryParse(str[0].Replace("RightLeg[", ""), System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positionx);
                    float.TryParse(str[1], System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positiony);
                    float.TryParse(str[2], System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positionz);
                    float.TryParse(str[3].Replace("]", ""), System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out scalarw);

                    System.Numerics.Quaternion quat = new System.Numerics.Quaternion();
                    quat.X = positionx;
                    quat.Y = positiony;
                    quat.Z = positionz;
                    quat.W = scalarw;


                    rightLegRotations[count] = quat;

                }
                else if (line.Contains("Neck"))
                {

                    float positionx, positiony, positionz, scalarw = 0.0f;
                    string[] str = line.Split(',');
                    float.TryParse(str[0].Replace("Neck[", ""), System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positionx);
                    float.TryParse(str[1], System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positiony);
                    float.TryParse(str[2], System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positionz);
                    float.TryParse(str[3].Replace("]", ""), System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out scalarw);

                    System.Numerics.Quaternion quat = new System.Numerics.Quaternion();
                    quat.X = positionx;
                    quat.Y = positiony;
                    quat.Z = positionz;
                    quat.W = scalarw;


                    neckRotations[count] = quat;

                }
                else if (line.Contains("LeftShoulder"))
                {

                    float positionx, positiony, positionz, scalarw = 0.0f;
                    string[] str = line.Split(',');
                    float.TryParse(str[0].Replace("LeftShoulder[", ""), System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positionx);
                    float.TryParse(str[1], System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positiony);
                    float.TryParse(str[2], System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positionz);
                    float.TryParse(str[3].Replace("]", ""), System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out scalarw);

                    System.Numerics.Quaternion quat = new System.Numerics.Quaternion();
                    quat.X = positionx;
                    quat.Y = positiony;
                    quat.Z = positionz;
                    quat.W = scalarw;


                    leftShoulderRotations[count] = quat;

                }
                else if (line.Contains("RightShoulder"))
                {

                    float positionx, positiony, positionz, scalarw = 0.0f;
                    string[] str = line.Split(',');
                    float.TryParse(str[0].Replace("RightShoulder[", ""), System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positionx);
                    float.TryParse(str[1], System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positiony);
                    float.TryParse(str[2], System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out positionz);
                    float.TryParse(str[3].Replace("]", ""), System.Globalization.NumberStyles.Any, CultureInfo.InvariantCulture, out scalarw);

                    System.Numerics.Quaternion quat = new System.Numerics.Quaternion();
                    quat.X = positionx;
                    quat.Y = positiony;
                    quat.Z = positionz;
                    quat.W = scalarw;


                    rightShoulderRotations[count] = quat;

                }

            }


        }

        readerStream.Close();


    /*   cursors[0] = new UnityEngine.Vector3(2, 1, -5);
        GameObject keyFrame = GameObject.CreatePrimitive(PrimitiveType.Sphere);
       keyFrame.transform.localScale = new UnityEngine.Vector3(0.1f, 0.1f, 0.1f);
       keyFrame.gameObject.name = "Cursor" + count;
       keyFrame.GetComponent<Renderer>().material.color = Color.red;

        keyFrame.gameObject.transform.position = cursors[0];
       
        cursors[1] = new UnityEngine.Vector3(6, 1, -5);
        GameObject keyFrame1 = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        keyFrame1.transform.localScale = new UnityEngine.Vector3(0.1f, 0.1f, 0.1f);
        keyFrame1.gameObject.name = "Cursor" + count;
        keyFrame1.GetComponent<Renderer>().material.color = Color.red;


        keyFrame1.gameObject.transform.position = cursors[1];

        count = 2;*/
       

    }

    public void Update()

        {
        
        //if cursor moves then the interpolation code runs.

        //1d code

     /*   if (gameObject.transform.hasChanged)
        {
            count = 0;
            Matrix<double> weights = calculateWeightCoefficients1D(cursors, rotations);

            double sum = calculateInterpolatedValue1D(0.0f, weights);            

            GameObject gm = GameObject.FindWithTag("LeftLeg");
            rotation.X = gm.transform.rotation.normalized.x;
            rotation.Y = gm.transform.rotation.normalized.y;
            rotation.Z = gm.transform.rotation.normalized.z;
            rotation.W = gm.transform.rotation.normalized.w;

            GameObject rightLeg = GameObject.FindWithTag("RightLeg");
            rotation.X = rightLeg.transform.rotation.normalized.x;
            rotation.Y = rightLeg.transform.rotation.normalized.y;
            rotation.Z = rightLeg.transform.rotation.normalized.z;
            rotation.W = rightLeg.transform.rotation.normalized.w;


            GameObject neck = GameObject.FindWithTag("Neck");
            rotation.X = neck.transform.rotation.normalized.x;
            rotation.Y = neck.transform.rotation.normalized.y;
            rotation.Z = neck.transform.rotation.normalized.z;
            rotation.W = neck.transform.rotation.normalized.w;


            Debug.Log("Cursor Position " + gameObject.transform.position.x);
            Debug.Log("sum "+ sum );

            //gm.transform.rotation = Quaternion.Euler((float)sum, 0, 0);
            //rightLeg.transform.rotation = Quaternion.Lerp(rightLeg.transform.rotation, Quaternion.Euler(0, 0, -sum), Time.deltaTime);
            //neck.transform.rotation = Quaternion.Lerp(rightLeg.transform.rotation, Quaternion.Euler(0, sum, 0), Time.deltaTime);




            gameObject.transform.hasChanged = false;

            count++;

        

        }*/

        //code for 3d rotations

        if(gameObject.transform.hasChanged)
        {

            GameObject leftLeg = GameObject.FindWithTag("LeftLeg");

            GameObject rightLeg = GameObject.FindWithTag("RightLeg");

            GameObject neck = GameObject.FindWithTag("Neck");

            GameObject leftShoulder = GameObject.FindWithTag("LeftShoulder");

            GameObject rightShoulder = GameObject.FindWithTag("RightShoulder");

         /*   System.Numerics.Quaternion firstPose = new System.Numerics.Quaternion();
            firstPose.X = 0.626317f;
            firstPose.Y = 0.1404756f;
            firstPose.Z = 0.7152255f;
            firstPose.W = 0.2764893f;


            System.Numerics.Quaternion secondPose = new System.Numerics.Quaternion();
            secondPose.X = 0.5487222f;
            secondPose.Y = 0.1102552f;
            secondPose.Z = 0.6743674f;
            secondPose.W = 0.4816397f;*/
      

            InterpolationFunction(rightShoulder, rightShoulderRotations);

            InterpolationFunction(neck, neckRotations);


            //rightLeg.transform.rotation = Quaternion.LerpUnclamped(rightLeg.transform.rotation.normalized, newPoseQuaternion.normalized, Time.deltaTime * 2);

            //for different joints the coefficient matrix needs to be inversed once.

            gameObject.transform.hasChanged = false;


        }

    }

    private void InterpolationFunction(GameObject rotObj, System.Numerics.Quaternion[] rotations)
    {
        //turn the quaternion into a rotation matrix.
        System.Numerics.Matrix4x4 rotationMatrixFirstPose = System.Numerics.Matrix4x4.CreateFromQuaternion(rotations[0]);
        System.Numerics.Matrix4x4 rotationMatrixSecondPose = System.Numerics.Matrix4x4.CreateFromQuaternion(rotations[1]);


        //System.Numerics.Matrix4x4 rotationMatrixThirdPose = System.Numerics.Matrix4x4.CreateFromQuaternion(rotations[2]);


        //find the interpolated entry for each entry of the matrix.
        System.Numerics.Matrix4x4 interpolatedMat = new System.Numerics.Matrix4x4();
        Vector3 position = gameObject.transform.position;


        interpolatedMat.M11 = interpolate(rotationMatrixFirstPose.M11, rotationMatrixSecondPose.M11, 0f, position);
        interpolatedMat.M12 = interpolate(rotationMatrixFirstPose.M12, rotationMatrixSecondPose.M12, 0f, position);
        interpolatedMat.M13 = interpolate(rotationMatrixFirstPose.M13, rotationMatrixSecondPose.M13, 0f, position);
        interpolatedMat.M21 = interpolate(rotationMatrixFirstPose.M21, rotationMatrixSecondPose.M21, 0f, position);
        interpolatedMat.M22 = interpolate(rotationMatrixFirstPose.M22, rotationMatrixSecondPose.M22, 0f, position);
        interpolatedMat.M23 = interpolate(rotationMatrixFirstPose.M23, rotationMatrixSecondPose.M23, 0f, position);
        interpolatedMat.M31 = interpolate(rotationMatrixFirstPose.M31, rotationMatrixSecondPose.M31, 0f, position);
        interpolatedMat.M32 = interpolate(rotationMatrixFirstPose.M32, rotationMatrixSecondPose.M32, 0f, position);
        interpolatedMat.M33 = interpolate(rotationMatrixFirstPose.M33, rotationMatrixSecondPose.M33, 0f, position);


        Debug.Log(interpolatedMat.ToString());

        Debug.Log("Cursor" + " " + position);

        //orthonormalize the interpolated Matrix

        System.Numerics.Matrix4x4 orthogonalMatrix = matrixOrthogonal(interpolatedMat);

        //Convert the matrix to quaternion.

        System.Numerics.Quaternion newPose = System.Numerics.Quaternion.CreateFromRotationMatrix(orthogonalMatrix);

        //apply the transformation to the joints

        UnityEngine.Quaternion newPoseQuaternion = new UnityEngine.Quaternion();
        newPoseQuaternion.x = newPose.X;
        newPoseQuaternion.y = newPose.Y;
        newPoseQuaternion.z = newPose.Z;
        newPoseQuaternion.w = newPose.W;


        rotObj.transform.rotation = Quaternion.Lerp(rotObj.transform.rotation, newPoseQuaternion, Time.deltaTime);
    }






    public System.Numerics.Matrix4x4 matrixOrthogonal(System.Numerics.Matrix4x4 matrix)
    {
     

     //Extract and orthogonalize vectors
        Vector3 v0 = new Vector3(matrix.M11, matrix.M12, matrix.M13);
        Vector3 v1 = new Vector3(matrix.M21, matrix.M22, matrix.M23);
        Vector3 v2 = new Vector3(matrix.M31, matrix.M32, matrix.M33);
        Vector3.OrthoNormalize(ref v0, ref v1, ref v2);


        System.Numerics.Matrix4x4 matrixOrtho = new System.Numerics.Matrix4x4();
        matrixOrtho.M11 = v0.x;
        matrixOrtho.M12 = v0.y;
        matrixOrtho.M13 = v0.z;
        matrixOrtho.M21 = v1.x;
        matrixOrtho.M22 = v1.y;
        matrixOrtho.M23 = v1.z;
        matrixOrtho.M31 = v2.x;
        matrixOrtho.M32 = v2.y;
        matrixOrtho.M33 = v2.z;
        matrixOrtho.M44 = 0.0F;

        // Debug.Log("Matrix: " + "\n" + matrix);
         Debug.Log("Matrix orthogonal: " + "\n" + matrixOrtho);

        return matrixOrtho;
    }




    private float interpolate(float m11, float m12,float m13, Vector3 position)
    {

        //this function has a limitation that it will work only for two points of interpolation.
        count = 2;
        
        if (count == 2)
        {
            double[,] delta = new double[6, 6];
            double[,] h = new double[6, 1];
            double[,] d = new double[6, 1];

            delta[0, 0] = 0;
            delta[0, 1] = findEulclideanDistance(cursors[0], cursors[1]);
            delta[1, 0] = findEulclideanDistance(cursors[0], cursors[1]);
            delta[1, 1] = 0;
            delta[1, 2] = 1;
            delta[0, 2] = 1;
            delta[2, 0] = 1;
            delta[2, 1] = 1;
            delta[0, 3] = cursors[0].x;
            delta[0, 4] = cursors[0].y;
            delta[0, 5] = cursors[0].z;
            delta[3, 0] = cursors[0].x;
            delta[4, 0] = cursors[0].y;
            delta[5, 0] = cursors[0].z;
            delta[1, 3] = cursors[1].x;
            delta[1, 4] = cursors[1].y;
            delta[1, 5] = cursors[1].z;
            delta[3, 1] = cursors[1].x;
            delta[4, 1] = cursors[1].y;
            delta[5, 1] = cursors[1].z;

            //  Debug.Log("deltaMatrix" + delta[0, 0] + " " + delta[0, 1] + " " + delta[0, 2] + " " + delta[0, 3] + " " + delta[1, 0] + " " + delta[1, 1] + " " + 
            // delta[1, 2] + " " + delta[1, 3] + " " + delta[2, 0] + " " + delta[2, 1] + " " + delta[2, 2] + " " + delta[2, 3] + " " + delta[3, 0] + " " + delta[3, 1] + " " + delta[3, 2] + " " + delta[3, 3]);     

            h[0, 0] = m11;
            h[1, 0] = m12;

            //  Debug.Log(" HMatrix" + h[0, 0] + " " + h[1, 0] + " "+ h[2,0] +" " + h[3,0] + " " + h[4,0]+ " "+ h[5,0]);

            //find the d and p0 p1 p2 and p3 for this entry by solving the linear system

            Matrix deltaMatrix = DenseMatrix.OfArray(delta);

        //    Debug.Log(deltaMatrix.ToString());

            Matrix hMat = DenseMatrix.OfArray(h);

         //   Debug.Log(hMat.ToString());

            Matrix<double> dMat = DenseMatrix.OfArray(d);

      
                deltaMatrix.TrySolveIterative(hMat, dMat, new MlkBiCgStab());
            

            Debug.Log(dMat.ToString());

           // Debug.Log(dMat[0, 0]);
         
            //getting the weights and the linear term now we have our interpolation function, we can get the final value for this entry by plugging in our gameObject or cursor position.

            double sum = dMat[0, 0] * findEulclideanDistance(position, cursors[0]) + dMat[1, 0] * findEulclideanDistance(position, cursors[1]) + dMat[2, 0] + dMat[3, 0] * position.x + dMat[4, 0] * position.y + dMat[5, 0] * position.z;
          
            return (float)sum;
         }

        //if three spatial keyframes are set.
        else if (count == 3)
        {

            double[,] delta = new double[7, 7];
            double[,] h = new double[7, 1];
            double[,] d = new double[7, 1];

            delta[0, 0] = 0;
            delta[0, 1] = findEulclideanDistance(cursors[0], cursors[1]);
            delta[0,2] = findEulclideanDistance(cursors[0], cursors[2]);
            delta[1, 0] = findEulclideanDistance(cursors[0], cursors[1]);
            delta[1, 1] = 0;
            delta[1, 2] = findEulclideanDistance(cursors[1], cursors[2]);
            delta[2, 0] = findEulclideanDistance(cursors[2], cursors[0]);
            delta[2, 1] = findEulclideanDistance(cursors[1], cursors[2]);
            delta[2, 2] = 0;


            delta[0, 3] = 1;
            delta[1, 3] = 1;
            delta[2, 3] = 1;
            delta[3, 0] = 1;
            delta[3, 1] = 1;
            delta[3, 2] = 1;


            delta[0, 4] = cursors[0].x;
            delta[0,5] =  cursors[0].y;
            delta[0,6] =  cursors[0].z;
            delta[1,4] = cursors[1].x;
            delta[1,5] = cursors[1].y;
            delta[1,6] = cursors[1].z;
            delta[2, 4] = cursors[2].x;
            delta[2, 5] = cursors[2].y;
            delta[2, 6] = cursors[2].z;


            delta[4,0] = cursors[0].x;
            delta[5,0] = cursors[0].y;
            delta[6,0] = cursors[0].z;
            delta[4,1] = cursors[1].x;
            delta[5,1] = cursors[1].y;
            delta[6,1] = cursors[1].z;
            delta[4,2] = cursors[2].x;
            delta[5,2] = cursors[2].y;
            delta[6, 2] = cursors[2].z;


            //  Debug.Log("deltaMatrix" + delta[0, 0] + " " + delta[0, 1] + " " + delta[0, 2] + " " + delta[0, 3] + " " + delta[1, 0] + " " + delta[1, 1] + " " + 
            // delta[1, 2] + " " + delta[1, 3] + " " + delta[2, 0] + " " + delta[2, 1] + " " + delta[2, 2] + " " + delta[2, 3] + " " + delta[3, 0] + " " + delta[3, 1] + " " + delta[3, 2] + " " + delta[3, 3]);     

            h[0, 0] = m11;
            h[1, 0] = m12;
            h[2, 0] = m13;

            //  Debug.Log(" HMatrix" + h[0, 0] + " " + h[1, 0] + " "+ h[2,0] +" " + h[3,0] + " " + h[4,0]+ " "+ h[5,0]);

            //find the d and p0 p1 p2 and p3 for this entry by solving the linear system

            Matrix deltaMatrix = DenseMatrix.OfArray(delta);

            Debug.Log(deltaMatrix.ToString());

            Matrix hMat = DenseMatrix.OfArray(h);

            Debug.Log(hMat.ToString());

            Matrix<double> dMat = DenseMatrix.OfArray(d);

            dMat = deltaMatrix.SolveIterative(hMat, new MlkBiCgStab());

            Debug.Log(dMat[0, 0]);

            // Debug.Log(dMat[0, 0]);

            //getting the weights and the linear term now we have our interpolation function, we can get the final value for this entry by plugging in our gameObject or cursor position.

            double sum = dMat[0, 0] * findEulclideanDistance(position, cursors[0]) + dMat[1, 0] * findEulclideanDistance(position, cursors[1]) + dMat[2,0] *findEulclideanDistance(position, cursors[2]) + dMat[3, 0] + dMat[4, 0] * position.x + dMat[5, 0] * position.y + dMat[6, 0] * position.z;

            return (float)sum;


        }





        throw new NotImplementedException();
    }

    private int findNearestCursorPosition(UnityEngine.Vector3[] cursors, UnityEngine.Vector3 position)
    {
        int index = 0;
        UnityEngine.Vector3 distanceVector = new UnityEngine.Vector3(0, 0, 0);
        for(int i = 0; i< count; i++)
        {
            UnityEngine.Vector3 temp = cursors[i] - position;

            if (temp.magnitude < distanceVector.magnitude)
            {
                index = i;
                distanceVector = temp;
            }
        }

        return index;
       
    }

    private double calculateInterpolatedValue1D(float matEntry, Matrix<double> weights)
    {

        double[,] result = new double[4, 1];
        Matrix<double> resultMat = DenseMatrix.OfArray(result);

        double sum = 0.0;

        for(int i =0; i< 2; i++)
        {

            sum += (float)weights[i, 0] * (findEulclideanDistance(gameObject.transform.position, cursors[i]));

        }

        if ((float)(weights[2, 0] + weights[3, 0] * gameObject.transform.position.x) != sum)   //po+ p1.x 
        {
            sum += (weights[2, 0] + weights[3, 0] * gameObject.transform.position.x);
        }

        return sum;
    }


    public float findEulclideanDistance(UnityEngine.Vector3 a, UnityEngine.Vector3 b)
    {

        return (float) Math.Abs(Math.Sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z)));
    }

    public Matrix<double> calculateWeightCoefficients1D(UnityEngine.Vector3[] cursors, System.Numerics.Quaternion[] rotations)
        {
        //function for the 1D examples in the final report

        //Convert the rotation into rotation matrix.
        int count = 2;

        //need a linear solver for the equation sum(delta(ci-cj))+polynomial function.


        double[,] h = new double[count + count, 1];
        double[,] d = new double[count + count, 1];
        double[,] result = new double[count + count, 1];
        double[,] delta = new double[count + count, count + count];


        //calculating the coefficient or the big matrix: need to calculate only once
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {

                if (i == count && j == count)
                    delta[i, j] = 0;
                else if (i == count)
                    delta[i, j] = 1;
                else if (j == count)
                    delta[i, j] = 1;
                else
                    delta[i, j] = Math.Abs(Math.Sqrt((cursors[i].x - cursors[j].x) * (cursors[i].x - cursors[j].x) + (cursors[i].y - cursors[j].y) * (cursors[i].y - cursors[j].y) + (cursors[i].z - cursors[j].z) * (cursors[i].z - cursors[j].z)));

            }

        }


        delta[0, 3] = cursors[0].x;
        delta[1, 3] = cursors[1].x;
        delta[3, 0] = cursors[0].x;
        delta[3, 1] = cursors[1].x;
        delta[1, 2] = 1;
        delta[2, 0] = 1;
        delta[2, 1] = 1;
        delta[2, 3] = 0;
        delta[3, 2] = 0;

        //can change the rotation values here to check result for different angles.
        h[0, 0] = (float)30.0f;

        h[1, 0] = (float)90.0f;

        //solver for the linear equation
        Matrix deltaMatrix = DenseMatrix.OfArray(delta);
        Matrix hMat = DenseMatrix.OfArray(h);
        Matrix dMat = DenseMatrix.OfArray(d);
        Matrix<double> resultMat = DenseMatrix.OfArray(result);

        Debug.Log("deltaMatrix" + delta[0, 0] + " " + delta[0, 1] + " " + delta[0, 2] + " " + delta[0, 3] + " " + delta[1, 0] + " " + delta[1, 1] + " " + delta[1, 2] + " " + delta[1, 3] + " " + delta[2, 0] + " " + delta[2, 1] + " " + delta[2, 2] + " " + delta[2, 3] + " " + delta[3, 0] + " " + delta[3, 1] + " " + delta[3, 2] + " " + delta[3, 3]);

        resultMat = deltaMatrix.Solve(hMat);

        Debug.Log(resultMat[0, 0] + " " + resultMat[1, 0] + " " + resultMat[2, 0] + " " + resultMat[3, 0]);

        return resultMat; //return the weighted coefficients found for these set of spatial keyframes



    }



}
