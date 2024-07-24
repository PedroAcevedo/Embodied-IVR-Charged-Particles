using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.Linq;


public class ParticleLines
{
    public GameObject[] particles;
    public float[] charges;
    public bool showForces;

    private float particleRadius = 0.25f;
    private Vector3[] lookUpTable = {
        new Vector3(0.0f, 5.0f, 0.0f),
        new Vector3(0.0f, 4.0f, 0.0f),
        new Vector3(0.0f, 3.0f, 0.0f),
        new Vector3(0.0f, 2.0f, 0.0f),
        new Vector3(0.0f, 1.0f, 0.0f),
        new Vector3(0.0f, 0.0f, 0.0f),
        new Vector3(0.0f, -1.0f, 0.0f),
        new Vector3(0.0f, -2.0f, 0.0f),
        new Vector3(0.0f, -3.0f, 0.0f),
        new Vector3(0.0f, -4.0f, 0.0f),
        new Vector3(0.0f, -5.0f, 0.0f),
    };
    private Vector2[] simulationLimits = {
        new Vector2(-9.5f, -9.5f),
        new Vector2(-9.5f, 9.5f),
        new Vector2(9.5f, 9.5f),
        new Vector2(9.5f, -9.5f),
    };

    private float lineDefaultWidth = 0.010f;
    private List<GameObject> lines = new List<GameObject>();
    private List<GameObject> arrows = new List<GameObject>();
    private int FIELD_LINES = 6;
    private float eps;
    private float EPSILON;
    private int linesLimit = 10;

    public ParticleLines(GameObject[] particles, float[] charges)
    {
        this.particles = particles;
        this.charges = charges;
        this.eps = particleRadius / 50;
        this.EPSILON = eps / 1.0E1f;
        this.showForces = true;
        Array.Sort(charges, particles);
    }

    
    public void Draw(bool mode)
    {
        CleanLines();

        for (int i = 0; i < particles.Length; i++)
        {
            if (mode)
            {
                if (isInsideSimulationBox2D(particles[i].transform.position))
                {
                    drawElectricLinesParticles2D(i);
                }
            }
            else
            {
                drawElectricLinesParticles3D(i);
            }
        }
    }

    public void CleanLines()
    {
        if (lines.Count > 0)
        {
            foreach (GameObject line in lines)
                GameObject.Destroy(line);
            foreach (GameObject arrow in arrows)
                GameObject.Destroy(arrow);
        }
    }

    public void AddArrow(Vector3 position, Vector3 nextPosition)
    {
        GameObject arrow = GameObject.Instantiate(Resources.Load("Prefabs/arrow"), position, Quaternion.identity) as GameObject;
        arrow.transform.LookAt(nextPosition);
        arrows.Add(arrow);

        Vector3 diff = position - nextPosition;
        diff.Normalize();

        float rot_z = Mathf.Atan2(diff.y, diff.x) * Mathf.Rad2Deg;
        arrow.transform.rotation = Quaternion.Euler(0f, 0f, rot_z - 180);
        arrow.transform.localScale -= new Vector3(0.0f, 0.00781458f, 0.0f);
    }

    private void drawElectricLinesParticles3D(int index)
    {
        for (int i = 0; i < FIELD_LINES; i++)
        {

            float x = particles[index].transform.position.x + eps * Mathf.Cos((float)(2 * Math.PI * i / FIELD_LINES));
            float y = particles[index].transform.position.y + eps * Mathf.Sin((float)(2 * Math.PI * i / FIELD_LINES));
            float z = particles[index].transform.position.z + eps * Mathf.Sin((float)(2 * Math.PI * i / FIELD_LINES));

            bool reachedAnotherCharge = false;

            // Check for infinite loop 
            bool infiniteLoop = false;
            int count = 0;
            float[] oldXs = { 0.0f, 0.0f };
            float[] oldYs = { 0.0f, 0.0f };
            float[] oldZs = { 0.0f, 0.0f };

            List<Vector3> lineField = new List<Vector3>();
            while (!reachedAnotherCharge && !infiniteLoop
                     && x > -linesLimit && x < linesLimit && y > -linesLimit && y < linesLimit
                     && z > -linesLimit && z < linesLimit)
            {

                // find the field (Ex, Ey, Ez) and field strength E at (x,y.z)
                float[] E = ETotal(x, y, z);
                float n = (float)Mathf.Sqrt(E[0] * E[0] + E[1] * E[1] + E[2] * E[2]);

                // if charge is negative the line needs to go backwards
                if (charges[index] > 0)
                {
                    x += E[0] / n * eps;
                    y += E[1] / n * eps;
                    z += E[2] / n * eps;
                }
                else
                {
                    x -= E[0] / n * eps;
                    y -= E[1] / n * eps;
                    z -= E[2] / n * eps;
                }

                lineField.Add(new Vector3(x, y, z));

                // stop in infinite loop
                if (Math.Abs(x - oldXs[0]) < EPSILON && Math.Abs(y - oldYs[0]) < EPSILON && Math.Abs(z - oldZs[0]) < EPSILON)
                {
                    infiniteLoop = true;
                }
                int index2 = count++ % 2;
                oldXs[index2] = x;
                oldYs[index2] = y;
                oldZs[index2] = z;


                // stop if the line ends in a charge
                for (int j = 0; j < charges.Length; j++)
                {
                    float dx = x - particles[j].transform.position.x;
                    float dy = y - particles[j].transform.position.y;
                    float dz = z - particles[j].transform.position.z;

                    if (Math.Sqrt(dx * dx + dy * dy + dz * dz) < eps) reachedAnotherCharge = true;
                }

            }

            AddNewLineRendererList(new Color(1.0f, 0.0f, 0.0f), lineField, charges[index]);

        }


    }

    private void drawElectricLinesParticles2D(int index)
    {
        for (int i = 0; i < FIELD_LINES; i++)
        {
            float x = particles[index].transform.position.x + eps * Mathf.Cos((float)(2 * Math.PI * i / FIELD_LINES));
            float y = particles[index].transform.position.y + eps * Mathf.Sin((float)(2 * Math.PI * i / FIELD_LINES));

            bool reachedAnotherCharge = false;

            // Check for infinite loop 
            bool infiniteLoop = false;
            int count = 0;
            float[] oldXs = { 0.0f, 0.0f };
            float[] oldYs = { 0.0f, 0.0f };

            List<Vector3> lineField = new List<Vector3>();
            while (!reachedAnotherCharge && !infiniteLoop
                     && x > -linesLimit && x < linesLimit && y > -linesLimit && y < linesLimit)
            {

                // find the field (Ex, Ey, Ez) and field strength E at (x,y.z)
                float[] E = ETotal2D(x, y);
                float n = (float)Mathf.Sqrt(E[0] * E[0] + E[1] * E[1]);

                // if charge is negative the line needs to go backwards
                if (charges[index] > 0)
                {
                    x += E[0] / n * eps;
                    y += E[1] / n * eps;
                }
                else
                {
                    x -= E[0] / n * eps;
                    y -= E[1] / n * eps;
                }

                lineField.Add(new Vector3(x, y, 0.3199998f));

                // stop in infinite loop
                if (Math.Abs(x - oldXs[0]) < EPSILON && Math.Abs(y - oldYs[0]) < EPSILON)
                {
                    infiniteLoop = true;
                }
                int index2 = count++ % 2;
                oldXs[index2] = x;
                oldYs[index2] = y;

                // stop if the line ends in a charge
                for (int j = 0; j < charges.Length; j++)
                {
                    float dx = x - particles[j].transform.position.x;
                    float dy = y - particles[j].transform.position.y;

                    if (Math.Sqrt(dx * dx + dy * dy) < eps) reachedAnotherCharge = true;
                }

            }

            AddNewLineRendererList(new Color(1.0f, 0.0f, 0.0f), lineField, charges[index]);

        }


    }

    private void AddNewLineRendererList(Color color, List<Vector3> pointList, float charge)
    {
        GameObject go = new GameObject($"LineRenderer_particle_1");
        LineRenderer goLineRenderer = go.AddComponent<LineRenderer>();
        goLineRenderer.material = new Material(Shader.Find("Sprites/Default"));
        goLineRenderer.startColor = Color.black;
        goLineRenderer.endColor = Color.black;
        goLineRenderer.startWidth = lineDefaultWidth;
        goLineRenderer.endWidth = lineDefaultWidth;

        goLineRenderer.positionCount = pointList.Count;
        goLineRenderer.SetPositions(pointList.ToArray());

        if (this.showForces)
        {
            int numberArrow = pointList.Count / 4;

            float distances = Vector2.Distance(new Vector2(pointList[0].x, pointList[0].y), new Vector2(pointList[pointList.Count - 1].x, pointList[pointList.Count - 1].y))
                + Vector2.Distance(new Vector2(pointList[0].x, pointList[0].y), new Vector2(pointList[(pointList.Count / 2)].x, pointList[(pointList.Count / 2)].y));

            if (inParticle(pointList[pointList.Count - 1], distances < 20))
            {
                numberArrow = pointList.Count / (20 - (int)distances);
            }

            if (charge > 0)
            {
                for (int i = 1; i < pointList.Count - numberArrow; i += numberArrow)
                {
                    AddArrow(pointList[i], pointList[i + 1]);
                }
            }
            else
            {
                for (int i = pointList.Count - numberArrow; i > 1; i -= numberArrow)
                {
                    AddArrow(pointList[i], pointList[i - 1]);
                }
            }

        }

        lines.Add(go);
    }

    private bool inParticle(Vector3 position, bool distances)
    {
        if (distances)
        {
            for (int i = 0; i < particles.Length; i++)
            {
                if (Vector2.Distance(new Vector2(particles[i].transform.position.x, particles[i].transform.position.y), new Vector2(position.x, position.y)) < 0.25)
                {
                    return true;
                }
            }
        }

        return false;
    }

    private float[] pointCharge(float charge, Vector3 position, float x, float y, float z)
    {
        float distance = (float)Math.Pow(Math.Pow(x - position.x, 2.0f) + Math.Pow(y - position.y, 2.0f) + Math.Pow(z - position.z, 2.0f), 0.5f);

        float[] chargeOnPoint = new float[3];

        chargeOnPoint[0] = charge * (x - position.x) / distance;
        chargeOnPoint[1] = charge * (y - position.y) / distance;
        chargeOnPoint[2] = charge * (z - position.z) / distance;

        return chargeOnPoint;
    }

    private float[] pointCharge2D(float charge, Vector3 position, float x, float y)
    {
        float distance = (float)Math.Pow(Math.Pow(x - position.x, 2.0f) + Math.Pow(y - position.y, 2.0f), 1.5f);

        float[] chargeOnPoint = new float[2];

        chargeOnPoint[0] = charge * (x - position.x) / distance;
        chargeOnPoint[1] = charge * (y - position.y) / distance;

        return chargeOnPoint;
    }

    private float[] ETotal(float x, float y, float z)
    {
        float[] Exy = new float[3];

        Exy[0] = 0.0f;
        Exy[1] = 0.0f;
        Exy[2] = 0.0f;

        for (int i = 0; i < charges.Length; i++)
        {
            float xp = particles[i].transform.position.x;
            float yp = particles[i].transform.position.y;
            float zp = particles[i].transform.position.z;

            if (xp > -linesLimit && xp < linesLimit && yp > -linesLimit && yp < linesLimit
                     && zp > -linesLimit && zp < linesLimit)
            {
                float[] E = pointCharge(charges[i], particles[i].transform.position, x, y, z);

                Exy[0] = Exy[0] + E[0];
                Exy[1] = Exy[1] + E[1];
                Exy[2] = Exy[2] + E[2];
            }
        }

        return Exy;
    }

    private float[] ETotal2D(float x, float y)
    {
        float[] Exy = new float[2];

        Exy[0] = 0.0f;
        Exy[1] = 0.0f;

        for (int i = 0; i < charges.Length; i++)
        {
            float[] E = pointCharge2D(charges[i], particles[i].transform.position, x, y);

            Exy[0] = Exy[0] + E[0];
            Exy[1] = Exy[1] + E[1];
        }

        return Exy;
    }

    private bool isInsideSimulationBox2D(Vector3 position)
    {
        Vector2 particlePos = new Vector2(position.x, position.y);

        float ABAM = Vector2.Dot(joinVector(simulationLimits[0], simulationLimits[1]), joinVector(simulationLimits[0], particlePos)); 
        float ABAB = Vector2.Dot(joinVector(simulationLimits[0], simulationLimits[1]), joinVector(simulationLimits[0], simulationLimits[1]));
        float BCBM = Vector2.Dot(joinVector(simulationLimits[1], simulationLimits[2]), joinVector(simulationLimits[1], particlePos));
        float BCBC = Vector2.Dot(joinVector(simulationLimits[1], simulationLimits[2]), joinVector(simulationLimits[1], simulationLimits[2]));

        return 0 <= ABAM && ABAM <= ABAB && 0 <= BCBM && BCBM <= BCBC;
    }

    private Vector2 joinVector(Vector2 p1, Vector2 p2)
    {
        return new Vector2(p2.x - p1.x, p2.y - p1.y);
    }
}
