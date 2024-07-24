using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class DataCollectionController
{

    private bool[] mappingButtons;

    public DataCollectionController()
    {
        mappingButtons = new bool[4];

        for (int i = 0; i < mappingButtons.Length; i++)
            mappingButtons[i] = false;
    }

    public void ButtonsPressed(ref SceneData scene)
    {
        if (OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger, OVRInput.Controller.RTouch) > 0.5 && !mappingButtons[0])
        {
            scene.rightHandButtonPress += 1;
            mappingButtons[0] = true;
        }
        else
        {
            if (OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger, OVRInput.Controller.RTouch) == 0 && mappingButtons[0])
            {
                mappingButtons[0] = false;
            }
        }

        if (OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger, OVRInput.Controller.RTouch) > 0.5 && !mappingButtons[1])
        {
            scene.rightHandButtonPress += 1;
            mappingButtons[1] = true;

        }
        else
        {
            if (OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger, OVRInput.Controller.RTouch) == 0 && mappingButtons[1])
            {
                mappingButtons[1] = false;
            }
        }

        // Left Hand

        if (OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger, OVRInput.Controller.LTouch) > 0.5 && !mappingButtons[2])
        {
            scene.leftHandButtonPress += 1;
            mappingButtons[2] = true;

        }
        else
        {
            if (OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger, OVRInput.Controller.LTouch) == 0 && mappingButtons[2])
            {
                mappingButtons[2] = false;
            }
        }

        if (OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger, OVRInput.Controller.LTouch) > 0.5 && !mappingButtons[3])
        {
            scene.leftHandButtonPress += 1;
            mappingButtons[3] = true;
        }
        else
        {
            if (OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger, OVRInput.Controller.LTouch) == 0 && mappingButtons[3])
            {
                mappingButtons[3] = false;
            }
        }


    }

    public void reportUser(ref SceneData scene, GameObject player)
    {
        scene.userPosition.Add(player.transform.position);
        scene.userRotation.Add(player.transform.rotation);
    }

}
