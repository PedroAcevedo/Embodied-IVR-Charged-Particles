using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class IndicatorController : MonoBehaviour
{

    public GameObject[] closeIndicators;
    public GameObject[] farIndicators;

    private bool showClose;

    // Start is called before the first frame update
    void Start()
    {
        showClose = true;
        IndicatorsActive(closeIndicators, true, false);
        IndicatorsActive(farIndicators, false, false);
    }

    // Update is called once per frame
    void Update()
    {
        if (showClose)
        {
            if (IndicatorsState(closeIndicators))
            {
                IndicatorsActive(closeIndicators, false, true);
                IndicatorsActive(farIndicators, true, true);
                showClose = false;
            }
        } 
        else
        {
            if (IndicatorsState(farIndicators))
            {
                IndicatorsActive(farIndicators, false, true);
            }
        }
    }


    bool IndicatorsState(GameObject[] indicators)
    {
        for (int i = 0; i < indicators.Length; ++i)
        {
            if(!indicators[i].GetComponent<IndicatorCollision>().isOccupied())
                return false;
        }
        return true;
    }

    void IndicatorReset(GameObject indicator)
    {
        indicator.GetComponent<IndicatorCollision>().removeParticleInteraction(false);
    }

    void IndicatorsActive(GameObject[] indicators, bool state, bool wait)
    {
        for (int i = 0; i < indicators.Length; ++i)
        {
            if (wait)
            {
                StartCoroutine(ExampleCoroutine(indicators[i], state));
            }
            else
            {
                indicators[i].SetActive(state);
            }
        }
    }

    IEnumerator ExampleCoroutine(GameObject indicator, bool state)
    {
        yield return new WaitForSeconds(2);
        indicator.SetActive(state);
        IndicatorReset(indicator);
    }

    void IndicatorsParticleReset(GameObject[] indicators)
    {
        for (int i = 0; i < indicators.Length; ++i)
        {
            indicators[i].GetComponent<IndicatorCollision>().cleanIndicator();
        }
    }

    public void resetIndicators()
    {
        showClose = true;
        IndicatorsActive(closeIndicators, true, false);
        IndicatorsParticleReset(closeIndicators);
        IndicatorsActive(farIndicators, false, false);
        IndicatorsParticleReset(farIndicators);
    }
}
