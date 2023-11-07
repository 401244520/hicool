


def train(model, train_loader, optimizer, criterion, device):
    model.train()
    train_loss = 0
    train_acc = 0
    for batch_idx, (data, target) in enumerate(train_loader):
        data, target = data.to(device), target.to(device)
        optimizer.zero_grad()
        output = model(data)
        loss = criterion(output, target)
        loss.backward()
        optimizer.step()
        train_loss += loss.item()
        pred = output.argmax(dim=1, keepdim=True)
        train_acc += pred.eq(target.view_as(pred)).sum().item()
    train_loss /= len(train_loader.dataset)
    train_acc /= len(train_loader.dataset)
    return train_loss, train_acc


def test(model, test_loader, criterion, device):
    model.eval()
    test_loss = 0
    test_acc = 0
    with torch.no_grad():
        for data, target in test_loader:
            data, target = data.to(device), target.to(device)
            output = model(data)
            test_loss += criterion(output, target).item()
            pred = output.argmax(dim=1, keepdim=True)
            test_acc += pred.eq(target.view_as(pred)).sum().item()
    test_loss /= len(test_loader.dataset)
    test_acc /= len(test_loader.dataset)
    return test_loss, test_acc

def plot_umap(embedding, labels, epoch):
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(embedding)
    plt.figure(figsize=(10, 10))
    plt.scatter(embedding[:, 0], embedding[:, 1], c=labels, cmap='jet')
    plt.colorbar()
    plt.title(f'UMAP visualization (epoch {epoch})')
    plt.show()
