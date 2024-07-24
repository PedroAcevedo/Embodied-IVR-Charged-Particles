using System.Collections;
using System.Collections.Generic;
using OVRTouchSample;
using UnityEngine;

public class CircleVibration : MonoBehaviour
{
    public GameObject RHand;            // VR right controller
    public GameObject LHand;            // VR left controller
    public float intensity;

    public float interactionTime = 0.0f;

    private float startTime = 0.0f;

    public void OnTriggerEnter(Collider other)
    { 
        if (other.gameObject.name == "RightHandAnchor")
        {
            OVRInput.SetControllerVibration(1, intensity, OVRInput.Controller.RTouch);
            startTime = Time.time;
        }

        if(other.gameObject.name == "LeftHandAnchor")
        {
            OVRInput.SetControllerVibration(1, intensity, OVRInput.Controller.LTouch);
            startTime = Time.time;
        }
    }
    public void OnTriggerExit(Collider other)
    {
        if (other.gameObject.name == "RightHandAnchor" || other.gameObject.name == "LeftHandAnchor")
        {
            interactionTime += (Time.time - startTime);
        }
    }
}
