using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using Unity.MathMatrix;

public class Sample : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {
        var a = new Vector2[] { new Vector2(5,5), new Vector2(15,5), new Vector2(15,15) , new Vector2(5,15) };
        var b = new Vector2[] { new Vector2(20,15), new Vector2(25,20), new Vector2(25,25), new Vector2(15,20) };

        var H = Matrix.DLT(a, b, false);
        H.DebugLog();

        var array = H.ToArray();
        Debug.Log(array[0]);
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
