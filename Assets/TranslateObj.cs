using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TranslateObj : MonoBehaviour
{
    
    private Vector3 screenPoint;
    private Vector3 offset;
    void OnMouseDown()
    {
        
            screenPoint = Camera.main.WorldToScreenPoint(gameObject.transform.position);

            offset = gameObject.transform.position - Camera.main.ScreenToWorldPoint(new Vector3(Input.mousePosition.x, Input.mousePosition.y, screenPoint.z));
            
    }

    [System.Obsolete]
    private void OnMouseDrag()
    {
 
        Vector3 curScreenPoint = new Vector3(Input.mousePosition.x, Input.mousePosition.y, screenPoint.z); // hardcode the y and z for your use

        Vector3 curPosition = Camera.main.ScreenToWorldPoint(curScreenPoint) + offset;
        gameObject.transform.position = curPosition;
       
    }
}

