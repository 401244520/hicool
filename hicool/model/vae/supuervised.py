import torch
import numpy as np
import umap
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

import torch
import torch.nn as nn
import torchvision.models as models

class ResNet(nn.Module):
    def __init__(self, input_size, num_classes):
        super(ResNet, self).__init__()
        self.input_size = input_size
        self.num_classes = num_classes
        self.resnet = models.resnet18(pretrained=True)
        self.resnet.conv1 = nn.Conv2d(1, 64, kernel_size=(7, 7), stride=(2, 2), padding=(3, 3), bias=False)
        self.resnet.fc = nn.Linear(512, num_classes)

    def forward(self, x):
        x = x.reshape(-1, 1, self.input_size, self.input_size)
        out = self.resnet(x)
        return out


# 定义训练函数
def train(model, data, labels, optimizer, criterion):
    model.train()
    optimizer.zero_grad()
    outputs = model(data)
    loss = criterion(outputs, labels)
    loss.backward()
    optimizer.step()
    return loss.item()

# 定义评估函数
def evaluate(model, data, labels):
    model.eval()
    outputs = model(data)
    preds = torch.argmax(outputs, dim=1)
    score = adjusted_rand_score(labels.numpy(), preds.numpy())
    return score

# 定义降维函数
def umap_embedding(data):
    reducer = umap.UMAP(random_state=42)
    embedding = reducer.fit_transform(data)
    return embedding

from 

# 定义超参数
input_size = data.shape[1]
output_size = len(torch.unique(labels))
num_epochs = 100
batch_size = 64
lr = 0.001

# 初始化模型、优化器和损失函数
model = ResNet(input_size, output_size)
optimizer = torch.optim.Adam(model.parameters(), lr=lr)
criterion = torch.nn.CrossEntropyLoss()

# 训练模型
for epoch in range(num_epochs):
    running_loss = 0.0
    for i in range(0, data.shape[0], batch_size):
        batch_data = data[i:i+batch_size]
        batch_labels = labels[i:i+batch_size]
        loss = train(model, batch_data, batch_labels, optimizer, criterion)
        running_loss += loss
    train_loss = running_loss / (data.shape[0] // batch_size)
    train_score = evaluate(model, data, labels)
    if epoch % 10 == 0:
        embedding = umap_embedding(data)
        plt.scatter(embedding[:, 0], embedding[:, 1], c=labels, cmap='Spectral')
        plt.show()
        print(f"Epoch {epoch}, Train Loss: {train_loss:.4f}, Train Score: {train_score:.4f}")

# 定义超参数
input_size = data.shape[1]
output_size = len(torch.unique(labels))
num_epochs = 100
batch_size = 64
lr = 0.001

# 初始化模型、优化器和损失函数
model = ResNet(input_size, output_size)
optimizer = torch.optim.Adam(model.parameters(), lr=lr)
criterion = torch.nn.CrossEntropyLoss()

# 定义训练函数
def train(model, dataloader, optimizer, criterion):
    model.train()
    running_loss = 0.0
    correct = 0
    total = 0
    for i, (inputs, labels) in enumerate(dataloader):
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
        running_loss += loss.item()
        _, predicted = torch.max(outputs.data, 1)
        total += labels.size(0)
        correct += (predicted == labels).sum().item()
    train_loss = running_loss / len(dataloader)
    train_acc = 100 * correct / total
    return train_loss, train_acc

# 定义可视化函数
def visualize(model, dataloader, epoch):
    model.eval()
    features = []
    labels = []
    with torch.no_grad():
        for inputs, target in dataloader:
            outputs = model(inputs)
            features.append(outputs.numpy())
            labels.append(target.numpy())
    features = np.concatenate(features, axis=0)
    labels = np.concatenate(labels, axis=0)
    umap_results = umap.UMAP(n_neighbors=5, min_dist=0.3, metric='euclidean').fit_transform(features)
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.scatter(umap_results[:, 0], umap_results[:, 1], c=labels, cmap='tab20')
    plt.title(f"UMAP visualization of epoch {epoch}")
    plt.savefig(f"UMAP_epoch{epoch}.png")
    plt.show()

# 训练和可视化
for epoch in range(num_epochs):
    train_loss, train_acc = train(model, trainloader, optimizer, criterion)
    if epoch % 10 == 9:
        visualize(model, trainloader, epoch+1)
    print(f"Epoch {epoch+1}/{num_epochs}: train_loss={train_loss:.4f}, train_acc={train_acc:.2f}%")
