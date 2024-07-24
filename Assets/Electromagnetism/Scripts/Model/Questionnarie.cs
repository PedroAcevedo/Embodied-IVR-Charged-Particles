using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[System.Serializable]
public class Questionnaire
{
    public List<SceneQuestions> scenes = new List<SceneQuestions>();
}

[System.Serializable]
public class SceneQuestions
{
    public string name;
    public List<Question> questions = new List<Question>();
}

[System.Serializable]
public class Question
{
    public List<string> options;
    public string body;
    public string figure;
    public string recordingFile = "";
    public bool hasImage;
    public string selectedAnswer = "";
}