using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class rotateObj : MonoBehaviour
{

    float rotateSpeed = 30;


    [System.Obsolete]
    private void OnMouseDrag()
    {
        float rotX = Input.GetAxis("Mouse X") * rotateSpeed * Mathf.Deg2Rad;
        float rotY = Input.GetAxis("Mouse Y") * rotateSpeed * Mathf.Deg2Rad;
        
        transform.RotateAround(Vector3.up, -rotX);
        transform.RotateAround(Vector3.right, rotY);
        
    }
}
