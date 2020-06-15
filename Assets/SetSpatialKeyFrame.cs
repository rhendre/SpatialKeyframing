using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;
using UnityEngine.Events;
using UnityEngine.EventSystems;

public class SetSpatialKeyFrame : MonoBehaviour

{
    int count = 0;

  
    public void OnSubmit()
    {
        //get the quaternian representation of the transfomation

        if (gameObject.name != "Cursor")
        {

           Quaternion rotation= gameObject.transform.rotation;
            
            //get the cursor of the representation

            Vector3 position = gameObject.transform.position;


            //find the rotation matrix from Quaternions.
            Debug.Log("OnSubmit called.");

            //save the spatial keyframe file which will further be used for interpolation

            string path =
              Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);

            // Write the string array to a new file named "WriteLines.txt".
            using (StreamWriter outputFile = File.AppendText(Path.Combine(path, "SpatialKeyFrame.txt")))
            {
                outputFile.WriteLine(gameObject.name +"["+ rotation.x + ", "+ rotation.y+", " + rotation.z + ", "+ rotation.w +"]");
                outputFile.Close();
            }
          


        }

        else
        {


            Vector3 position = gameObject.transform.position;


            GameObject keyFrame = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            keyFrame.transform.localScale = new Vector3(0.1f, 0.1f, 0.1f);
            keyFrame.gameObject.name = "Cursor" + count;
            count++;

            keyFrame.gameObject.transform.position = position;


            //find the rotation matrix from Quaternions.
            Debug.Log("OnSubmit called.");

            //save the spatial keyframe file which will further be used for interpolation

            string path =
              Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);

            // Write the string array to a new file named "WriteLines.txt".
            using (StreamWriter outputFile = File.AppendText(Path.Combine(path, "SpatialKeyFrame.txt")))
            {
                outputFile.WriteLine("Cursor[" + position.x + ", " + position.y + ", " + position.z + "]");
                outputFile.Close();
            }


            count++;

        }


        




      }



    

}


