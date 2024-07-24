using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System;

public class UserReportController
{
    private List<SceneData> scenes = new List<SceneData>();
    private List<UserResponse> userResponses = new List<UserResponse>();

    public void SceneInfo(SceneData currentScene)
    {
        scenes.Add(currentScene);
    }

    public void generateUserResponses()
    {
        foreach (SceneQuestions scene in QuizController.questionnaire.scenes)
        {
            List<string> responses = new List<string>();

            for (int i = 0;  i < scene.questions.Count; i++)
            {
                if (scene.questions[i].recordingFile != "")
                {
                    responses.Add(scene.questions[i].recordingFile);
                }
                else
                {
                    responses.Add(scene.questions[i].selectedAnswer + "");
                }
            }

            userResponses.Add(new UserResponse(responses, scene.name));
        }
    }

    public void SaveIntoJson(string userId, string recording)
    {
        ScenesData levelsData = new ScenesData(userId, !SceneController.LeftHander, recording, scenes, userResponses);
        string data = JsonUtility.ToJson(levelsData);
        System.IO.File.WriteAllText(Application.persistentDataPath + "/user-" + userId + "-" + DateTimeOffset.Now.ToUnixTimeMilliseconds() + "-.json", data);
    }

}