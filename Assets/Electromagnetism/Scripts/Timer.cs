using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Timer
{
    public float currentTime;
    public bool timerIsRunning = false;

    public Timer()
    {
        this.currentTime = 0;
    }

    public void start()
    {
        this.timerIsRunning = true;
    }

    public void run()
    {
        if (timerIsRunning)
        {
            currentTime += Time.deltaTime;
        }
    }

    public string DisplayTimeInMinutes()
    {
        float minutes = Mathf.FloorToInt(currentTime / 60);
        float seconds = Mathf.FloorToInt(currentTime % 60);

        return string.Format("{0:00}:{1:00}", minutes, seconds);
    }


    public string DisplayTimeInSeconds()
    {
        return string.Format("{0:00}", currentTime);
    }

    public void reset()
    {
        this.currentTime = 0;
        this.timerIsRunning = true;
    }

    public void stop()
    {
        this.timerIsRunning = false;
    }
}
