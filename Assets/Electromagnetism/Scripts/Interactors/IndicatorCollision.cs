using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class IndicatorCollision : MonoBehaviour
{
    public Vector3 particlePosition;
    public GameObject currentParticle;

    public void OnTriggerStay(Collider other)
    {
        if (other.gameObject.name.Contains("Particle"))
        {
            currentParticle = other.gameObject;
            removeParticleInteraction(true);
            currentParticle.transform.position = particlePosition;
        }
    }

    public void removeParticleInteraction(bool value)
    {
        if(currentParticle != null)
        {
            currentParticle.GetComponent<Rigidbody>().freezeRotation = value;

            if (value)
            {
                currentParticle.GetComponent<Rigidbody>().constraints = RigidbodyConstraints.FreezePosition;
            }
            else
            {
                currentParticle.GetComponent<Rigidbody>().constraints = RigidbodyConstraints.None;
                currentParticle.GetComponent<Rigidbody>().constraints = RigidbodyConstraints.FreezePositionZ | RigidbodyConstraints.FreezeRotationZ;
            }
        }
    }

    void OnTriggerExit(Collider other)
    {
        if (currentParticle != null)
        {
            currentParticle.transform.position = particlePosition;
        }
    }

    public bool isOccupied()
    {
        if (currentParticle != null)
        {
            if (!currentParticle.GetComponent<OVRGrabbable>().isGrabbed)
            {
                return true;
            }
        }

        return false;
    }

    public void cleanIndicator()
    {
        currentParticle = null;
    }
}
