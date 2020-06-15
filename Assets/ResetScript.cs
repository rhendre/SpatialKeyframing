using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ResetScript : MonoBehaviour
{
    public void OnReset()
    {
        if (gameObject.name == "Neck")
        {
            Quaternion neck = new Quaternion(-0.6f, 0.5f, -0.5f, 0.4f);
            gameObject.transform.rotation = neck;
        }
         if(gameObject.name == "LeftLeg")
        {    
            Quaternion rightLeg = new Quaternion(0.6f, 0.4f, 0.6f, 0.4f);
            gameObject.transform.rotation = rightLeg;
        }
         if(gameObject.name == "RightLeg")

        {
            Quaternion rightLeg = new Quaternion(-0.6f, 0.5f, -0.5f, 0.4f);
            gameObject.transform.rotation = rightLeg;
        }
    }
}
