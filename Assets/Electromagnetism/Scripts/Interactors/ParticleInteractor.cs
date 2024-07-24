using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Events;
using UnityEngine.EventSystems;

public class ParticleInteractor : MonoBehaviour
{
    private RaycastHit hit;
    private GameObject target;

    // Update is called once per frame
    void Update()
    {
        if (Physics.Raycast(this.transform.position, this.transform.forward, out hit))
        {
            if (hit.transform.gameObject.name.Contains("Particle") )
            {
                if (target != hit.collider.gameObject) 
                {
                    if (target != null) 
                    {
                        target.transform.localScale = new Vector3(0.5f, 0.5f, 0.5f);
                    }
                    target = hit.collider.gameObject; // make object my target
                    target.transform.localScale = new Vector3(0.7f, 0.7f, 0.7f);
                    Debug.Log("Change Size");
                }
            } else
            {
                if (target != null)
                {
                    target.transform.localScale = new Vector3(1f, 1f, 1f);
                }
                target = null;
            }
        }
    }

}
