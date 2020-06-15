using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;
using UnityEngine.Events;
using UnityEngine.EventSystems;

public class Script : EventTrigger
{


    public override void OnSubmit(BaseEventData data)
    {
        //get the quaternian representation of the transfomation

        UnityEngine.Quaternion rot = gameObject.transform.rotation;

        //get the cursor of the representation

        Vector3 cursor = gameObject.transform.position;

        string name = gameObject.name;


        //find the rotation matrix from Quaternions.
        Debug.Log("OnSubmit called.");

        //save the spatial keyframe file which will further be used for interpolation

        string path =
          Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);

        // Write the string array to a new file named "WriteLines.txt".
        using (StreamWriter outputFile = File.AppendText(Path.Combine(path, "SpatialKeyFrame.txt")))
        {
            outputFile.WriteLine(name + rot + cursor);
        }


    }

}


